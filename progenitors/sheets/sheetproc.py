"""Google Sheets API integration for progenitor data download/upload.

Settings: sheet column names and OAuth scopes in progenitors/settings/sheets.py
"""
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from astropy.table import Table
import io
import logging
import os
import copy
import pickle
import random
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import zstandard as zstd

from ..settings.sheets import SHEET_COLUMNS, SCOPES

logger = logging.getLogger(__name__)

# Per-sheet Table.meta cache: zstd-compressed pickle (smaller/faster on slow FS).
METADATA_FILE_SUFFIX = ".pkl.zst"
LEGACY_METADATA_SUFFIX = ".pkl"
METADATA_BACKUP_SUFFIX = ".backup.pkl.zst"
LEGACY_BACKUP_SUFFIX = ".backup.pkl"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"

# Default Sheets API v4 write quota is 60/min per user; pace bursts + retry 429/503.
_DEFAULT_MIN_WRITE_INTERVAL = float(os.environ.get("PROGENITORS_SHEETS_WRITE_INTERVAL", "1.05"))
_last_write_time = 0.0


def _pace_sheets_writes():
    """Space write calls to reduce 429 WriteRequestsPerMinutePerUser bursts."""
    global _last_write_time
    now = time.monotonic()
    wait = _DEFAULT_MIN_WRITE_INTERVAL - (now - _last_write_time)
    if wait > 0:
        time.sleep(wait)
    _last_write_time = time.monotonic()


def _http_status(err):
    r = getattr(err, "resp", None)
    if r is None:
        return None
    return getattr(r, "status", None)


def _retry_after_seconds(err):
    r = getattr(err, "resp", None)
    if r is None or not hasattr(r, "get"):
        return None
    for key in ("retry-after", "Retry-After"):
        val = r.get(key)
        if val is not None:
            try:
                return float(val)
            except (TypeError, ValueError):
                continue
    return None


def _execute_sheets_read(request_fn, max_attempts=8, base_delay=1.0):
    """Run a Sheets read (get, values.get); retry on 429/503 without write pacing."""
    for attempt in range(max_attempts):
        try:
            return request_fn()
        except HttpError as e:
            status = _http_status(e)
            if status not in (429, 503):
                raise
            delay = base_delay * (2 ** attempt) + random.uniform(0, 0.5)
            ra = _retry_after_seconds(e)
            if ra is not None:
                delay = max(delay, ra)
            logger.warning(
                "Google Sheets read %s; sleeping %.1fs (retry %d/%d)",
                status,
                delay,
                attempt + 1,
                max_attempts,
            )
            time.sleep(delay)
    raise RuntimeError("Google Sheets API: read retries exhausted")


def _execute_sheets_write(request_fn, max_attempts=12, base_delay=2.0):
    """
    Run a Sheets write (batchUpdate, values.update, etc.); retry on rate limit.

    Handles HTTP 429 and 503 with exponential backoff and optional Retry-After.
    """
    for attempt in range(max_attempts):
        _pace_sheets_writes()
        try:
            return request_fn()
        except HttpError as e:
            status = _http_status(e)
            if status not in (429, 503):
                raise
            delay = base_delay * (2 ** attempt) + random.uniform(0, 1.0)
            ra = _retry_after_seconds(e)
            if ra is not None:
                delay = max(delay, ra)
            logger.warning(
                "Google Sheets write %s (rate limit or unavailable); sleeping %.1fs then retry %d/%d",
                status,
                delay,
                attempt + 1,
                max_attempts,
            )
            time.sleep(delay)
    raise RuntimeError("Google Sheets API: exceeded retries after rate limits")

# Repo root: from progenitors/sheets/sheetproc.py go up to repo
_here = os.path.abspath(os.path.dirname(__file__))
basedir = os.path.dirname(os.path.dirname(_here))
_credentials = os.environ.get("PROGENITORS_CREDENTIALS_PATH") or os.path.join(basedir, "credentials.json")

params = {
    "SHEET": os.environ.get("PROGENITORS_SHEET", ""),
    "token": os.path.join(basedir, "token.pickle"),
    "credentials": _credentials,
    "target": basedir,
    "yse": {
        "user": os.environ.get("YSE_USER", ""),
        "password": os.environ.get("YSE_PASSWORD", ""),
    },
    "cols": list(SHEET_COLUMNS),
}


def transpose(data):
    return list(map(list, zip(*data)))


def _atomic_pickle_dump(obj, path, protocol=None):
    """
    Pickle ``obj`` to ``path`` atomically: write a temp file in the same
    directory, fsync, then ``os.replace`` into place so an interrupted or
    failed write never truncates an existing ``path``.
    """
    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL
    path = os.path.abspath(path)
    directory = os.path.dirname(path) or os.getcwd()
    os.makedirs(directory, exist_ok=True)
    basename = os.path.basename(path)
    fd, tmp_path = tempfile.mkstemp(
        prefix=f".{basename}.", suffix=".tmp", dir=directory
    )
    try:
        with os.fdopen(fd, "wb") as f:
            pickle.dump(obj, f, protocol=protocol)
            f.flush()
            try:
                os.fsync(f.fileno())
            except OSError:
                pass
        os.replace(tmp_path, path)
    except BaseException:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        except OSError:
            pass
        raise


def _metadata_zstd_level() -> int:
    """Default 1 (fast); raise for smaller files (e.g. 3–6) via env."""
    try:
        return int(os.environ.get("PROGENITORS_METADATA_ZSTD_LEVEL", "1"))
    except ValueError:
        return 1


def _metadata_decompress_max_bytes() -> int:
    """Upper bound for one-shot zstd output; avoids 'Destination buffer is too small'."""
    default = 2 * 1024**3
    try:
        n = int(os.environ.get("PROGENITORS_METADATA_MAX_DECOMPRESS_BYTES", str(default)))
    except ValueError:
        n = default
    return max(n, 1)


# zstd contexts are not thread-safe; load/save metadata uses ThreadPoolExecutor.
_zstd_tls = threading.local()


def _thread_zstd_compressor(level: int) -> zstd.ZstdCompressor:
    d = getattr(_zstd_tls, "compressors", None)
    if d is None:
        d = {}
        _zstd_tls.compressors = d
    if level not in d:
        d[level] = zstd.ZstdCompressor(level=level)
    return d[level]


def _thread_zstd_decompressor() -> zstd.ZstdDecompressor:
    d = getattr(_zstd_tls, "decompressor", None)
    if d is None:
        d = zstd.ZstdDecompressor()
        _zstd_tls.decompressor = d
    return d


def _compress_metadata_obj(obj: object) -> bytes:
    level = _metadata_zstd_level()
    raw = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    return _thread_zstd_compressor(level).compress(raw)


def _decompress_zstd_to_bytes(dctx: zstd.ZstdDecompressor, blob: bytes) -> bytes:
    """Decompress zstd; large pickles need max_output_size (default one-shot buffer is small)."""
    max_out = _metadata_decompress_max_bytes()
    last_zerr: zstd.ZstdError | None = None
    for kwargs in ({"max_output_size": max_out}, {}):
        try:
            return dctx.decompress(blob, **kwargs)
        except TypeError:
            continue
        except zstd.ZstdError as e:
            last_zerr = e
            break
    if last_zerr is not None:
        low = str(last_zerr).lower()
        if "too small" in low or "destination buffer" in low:
            reader = dctx.stream_reader(io.BytesIO(blob))
            out = bytearray()
            chunk = 8 * 1024 * 1024
            while True:
                block = reader.read(chunk)
                if not block:
                    break
                out.extend(block)
            return bytes(out)
        raise last_zerr
    raise RuntimeError("zstd decompress: no compatible decompress() signature")


def _decompress_metadata_blob(blob: bytes) -> object:
    if len(blob) >= 4 and blob[:4] == _ZSTD_MAGIC:
        raw = _decompress_zstd_to_bytes(_thread_zstd_decompressor(), blob)
    else:
        raw = blob
    return pickle.loads(raw)


def _atomic_metadata_zstd_dump(obj: object, path: str) -> None:
    """Atomic write of zstd-compressed pickle (same replace semantics as _atomic_pickle_dump)."""
    blob = _compress_metadata_obj(obj)
    path = os.path.abspath(path)
    directory = os.path.dirname(path) or os.getcwd()
    os.makedirs(directory, exist_ok=True)
    basename = os.path.basename(path)
    fd, tmp_path = tempfile.mkstemp(
        prefix=f".{basename}.", suffix=".tmp", dir=directory
    )
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(blob)
            f.flush()
            if os.environ.get("PROGENITORS_METADATA_FSYNC", "1").strip().lower() not in (
                "0",
                "false",
                "no",
            ):
                try:
                    os.fsync(f.fileno())
                except OSError:
                    pass
        os.replace(tmp_path, path)
    except BaseException:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        except OSError:
            pass
        raise


def _metadata_savekey(key):
    return key.replace(" ", "_").replace("-", "_").replace("/", "_")


def _load_one_metadata_file(key, path):
    try:
        with open(path, "rb") as f:
            blob = f.read()
        meta = _decompress_metadata_blob(blob)
        return key, meta, None
    except (
        EOFError,
        MemoryError,
        pickle.UnpicklingError,
        OSError,
        zstd.ZstdError,
        ValueError,
    ) as e:
        return key, None, (path, e)


def _metadata_fallback_chain(directory: str, savekey: str) -> list[tuple[str, str]]:
    """Ordered (path, label) candidates; only existing paths are tried by the loader."""
    return [
        (os.path.join(directory, savekey + METADATA_FILE_SUFFIX), "main"),
        (
            os.path.join(directory, "backup", savekey + METADATA_BACKUP_SUFFIX),
            "backup",
        ),
        (os.path.join(directory, savekey + LEGACY_METADATA_SUFFIX), "legacy"),
        (
            os.path.join(directory, "backup", savekey + LEGACY_BACKUP_SUFFIX),
            "legacy_backup",
        ),
    ]


def _metadata_any_candidate_exists(directory: str, savekey: str) -> bool:
    return any(
        os.path.isfile(p) for p, _ in _metadata_fallback_chain(directory, savekey)
    )


def _maybe_quarantine_corrupt_primary(primary_path: str) -> None:
    """
    After a successful load from backup/legacy, move an unreadable primary .pkl.zst
    aside so the next run does not retry decompression on it every time.
    Disable with PROGENITORS_METADATA_QUARANTINE_CORRUPT=0.
    """
    if not os.path.isfile(primary_path):
        return
    if os.environ.get("PROGENITORS_METADATA_QUARANTINE_CORRUPT", "1").strip().lower() in (
        "0",
        "false",
        "no",
    ):
        return
    bad = primary_path + ".corrupt"
    try:
        if os.path.isfile(bad):
            os.remove(bad)
        os.replace(primary_path, bad)
        logger.warning(
            "metadata: quarantined corrupt primary to %s (loaded from backup/legacy)",
            os.path.basename(bad),
        )
    except OSError as e:
        logger.warning("metadata: could not quarantine %s: %s", primary_path, e)


def _load_metadata_with_fallbacks(key, directory: str, savekey: str):
    """
    Try main .pkl.zst, then backup .backup.pkl.zst, then legacy .pkl, then
    legacy backup. Returns (key, meta, None) on success, or (key, None, errors).
    """
    primary = os.path.join(directory, savekey + METADATA_FILE_SUFFIX)
    errors = []
    for path, label in _metadata_fallback_chain(directory, savekey):
        if not os.path.isfile(path):
            continue
        _k, meta, err = _load_one_metadata_file(key, path)
        if err is None:
            if label != "main":
                logger.warning(
                    "metadata: recovered %r from %s (%s); primary cache was unreadable. "
                    "Run the pipeline once to rewrite %s. Common causes: cloud sync "
                    "(e.g. Dropbox) serving a partial file, interrupted save, or "
                    "conflicted copy — check the host app’s version history if needed.",
                    savekey,
                    os.path.basename(path),
                    label,
                    METADATA_FILE_SUFFIX,
                )
                if os.path.isfile(primary) and path != primary:
                    _maybe_quarantine_corrupt_primary(primary)
            return (key, meta, None)
        errors.append(err)
    return (key, None, errors)


def initiate_gsheet(token_file, credentials=None):
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as token:
            creds = pickle.load(token)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(credentials, SCOPES)
            creds = flow.run_local_server(port=0)
        _atomic_pickle_dump(creds, token_file)
    service = build('sheets', 'v4', credentials=creds)
    return service.spreadsheets()


def download_progenitor_data(spreadsheetId):
    sheet = initiate_gsheet(params['token'], credentials=params['credentials'])
    response = sheet.get(spreadsheetId=spreadsheetId).execute()
    if not response or 'sheets' not in response:
        return {}
    alldata = {}
    for subsheet in response['sheets']:
        name = subsheet['properties']['title']
        data = sheet.values().get(spreadsheetId=spreadsheetId, range=name).execute()
        if 'values' in data:
            header = data['values'][0]
            bodydata = data['values'][1:]
            bodydata = [b for b in bodydata if len(b[0].strip()) > 0]
            tabdata = transpose(bodydata)
            if len(tabdata) == 0:
                tabdata = [['']] * len(header)
            elif len(tabdata) < len(header):
                for _ in range(len(header) - len(tabdata)):
                    tabdata.append([' ' * 100] * len(tabdata[0]))
            table = Table(tabdata, names=header)
            table.meta['name'] = name
            remove_rows = [i for i, row in enumerate(table) if len(row['Name'].strip()) == 0]
            table.remove_rows(remove_rows)
            alldata[name] = table
        else:
            init_data = [['X' * 100] for _ in params['cols']]
            table = Table(init_data, names=params['cols'])
            alldata[name] = table
    return alldata


def convert_table_to_lists(table, mask=None):
    output = [list(params['cols'])]
    for i, row in enumerate(table):
        if mask is not None and mask[i]:
            continue
        add_row = []
        for k in output[0]:
            if k in row.colnames:
                add_row.append('\'' + str(row[k]))
            elif k.lower() in row.colnames:
                add_row.append('\'' + str(row[k.lower()]))
            else:
                add_row.append('')
        output.append(add_row)
    return output


def upload_progenitor_data(spreadsheetId, all_tables, mask=False):
    """
    Upload table(s) to the spreadsheet.

    Pass only the tab(s) you changed when possible: each tab triggers several
    API writes; uploading the full workbook every iteration exceeds the default
    60 writes/minute/user quota.
    """
    sheet = initiate_gsheet(params['token'])
    response = _execute_sheets_read(
        lambda: sheet.get(spreadsheetId=spreadsheetId).execute()
    )
    if not response or 'sheets' not in response:
        return
    for key in all_tables.keys():
        mask_vals = all_tables[key].meta.get('mask') if mask else None
        table_data = convert_table_to_lists(all_tables[key], mask=mask_vals)
        sheetid_list = [s['properties']['sheetId'] for s in response['sheets']
                       if s['properties']['title'] == key]
        if len(sheetid_list) == 0:
            body = {'requests': [{'addSheet': {'properties': {'title': key}}}]}
            _execute_sheets_write(
                lambda b=body: sheet.batchUpdate(
                    spreadsheetId=spreadsheetId, body=b
                ).execute()
            )
            response = _execute_sheets_read(
                lambda: sheet.get(spreadsheetId=spreadsheetId).execute()
            )
            sheetid = [s['properties']['sheetId'] for s in response['sheets']
                       if s['properties']['title'] == key][0]
        else:
            sheetid = sheetid_list[0]
        # One batchUpdate instead of two (fewer write quota units per tab).
        prep_body = {
            'requests': [
                {
                    'updateDimensionProperties': {
                        'range': {
                            'sheetId': sheetid,
                            'dimension': 'ROWS',
                            'startIndex': 0,
                            'endIndex': 9999,
                        },
                        'properties': {'hiddenByUser': False},
                        'fields': 'hiddenByUser',
                    }
                },
                {
                    'updateCells': {
                        'range': {'sheetId': sheetid},
                        'fields': 'userEnteredValue',
                    }
                },
            ]
        }
        _execute_sheets_write(
            lambda b=prep_body: sheet.batchUpdate(
                spreadsheetId=spreadsheetId, body=b
            ).execute()
        )
        values_body = {'majorDimension': 'ROWS', 'values': table_data}
        _execute_sheets_write(
            lambda vb=values_body: sheet.values().update(
                spreadsheetId=spreadsheetId,
                range=key,
                valueInputOption='USER_ENTERED',
                body=vb,
            ).execute()
        )


def load_metadata_from_file(all_data, directory):
    loaded = []
    errors = []
    items = []
    for key in all_data.keys():
        savekey = _metadata_savekey(key)
        if _metadata_any_candidate_exists(directory, savekey):
            items.append((key, savekey, directory))

    if not items:
        logger.info(
            "metadata load: no %s / backup / legacy in %s (starting fresh meta)",
            METADATA_FILE_SUFFIX,
            directory,
        )
        return all_data

    if len(items) == 1:
        results = [
            _load_metadata_with_fallbacks(items[0][0], items[0][2], items[0][1])
        ]
    else:
        max_workers = min(16, len(items))
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            results = list(
                ex.map(
                    lambda t: _load_metadata_with_fallbacks(t[0], t[2], t[1]),
                    items,
                )
            )

    for (key, savekey, _directory), (_k, meta, err_list) in zip(items, results):
        if meta is not None:
            all_data[key].meta = meta
            loaded.append(savekey)
        else:
            for path, err in err_list or []:
                errors.append((path, err))

    if loaded:
        logger.info(
            "metadata load: %d sheet(s) from %s — %s",
            len(loaded),
            directory,
            ", ".join(f"{sk}.pkl.zst" for sk in loaded),
        )
    for path, err in errors:
        logger.warning("metadata load failed %s: %s", path, err)
    return all_data


def save_metadata_to_file(all_data, directory, make_backup=False):
    os.makedirs(directory, exist_ok=True)
    if make_backup:
        os.makedirs(os.path.join(directory, "backup"), exist_ok=True)

    rows = []
    for key in all_data.keys():
        savekey = _metadata_savekey(key)
        fulloutname = os.path.join(directory, savekey + METADATA_FILE_SUFFIX)
        backupname = os.path.join(directory, "backup", savekey + METADATA_BACKUP_SUFFIX)
        legacy_main = os.path.join(directory, savekey + LEGACY_METADATA_SUFFIX)
        legacy_backup = os.path.join(directory, "backup", savekey + LEGACY_BACKUP_SUFFIX)
        rows.append((key, fulloutname, backupname, legacy_main, legacy_backup))

    def _save_row(row):
        key, fulloutname, backupname, legacy_main, legacy_backup = row
        _atomic_metadata_zstd_dump(all_data[key].meta, fulloutname)
        if make_backup:
            _atomic_metadata_zstd_dump(all_data[key].meta, backupname)
        for leg in (legacy_main, legacy_backup):
            if os.path.isfile(leg):
                try:
                    os.remove(leg)
                except OSError:
                    logger.debug("metadata: could not remove legacy file %s", leg)
        return fulloutname

    if len(rows) <= 1:
        paths = [_save_row(r) for r in rows]
    else:
        max_workers = min(16, len(rows))
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            paths = list(ex.map(_save_row, rows))

    if paths:
        logger.info(
            "metadata save: %d compressed cache file(s)%s — %s",
            len(paths),
            " +backup" if make_backup else "",
            "; ".join(os.path.basename(p) for p in paths),
        )
