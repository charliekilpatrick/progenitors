"""
Main pipeline: build progenitor data from Google Sheets and metadata sources.

Downloads sheet data, fetches metadata (NED, TNS, YSE, ADS, OSC, MAST, etc.),
optionally updates classification and sends email alerts for newly binned core
types. Settings (default
metadata list, trim keys, YSE query) are in progenitors.settings.pipeline.
"""
import logging
import os
import time
import copy
from astropy.table import unique, Table

from . import util
from .sheets import sheetproc
from .settings.pipeline import (
    TRIM_METADATA_KEEP_KEYS,
    DEFAULT_METADATA,
    DEFAULT_YSE_SQL_QUERY,
    ALERT_RECLASSIFICATION_SHEET_KEYS,
)

logger = logging.getLogger(__name__)


class _Profiler:
    def __init__(self, enabled: bool):
        self.enabled = enabled
        self._marks = {}
        self.records = []

    def start(self, key: str) -> None:
        if self.enabled:
            self._marks[key] = time.perf_counter()

    def stop(self, key: str, label: str) -> None:
        if self.enabled and key in self._marks:
            dt = time.perf_counter() - self._marks.pop(key)
            self.records.append((label, dt))

    def log_summary(self) -> None:
        if not self.enabled:
            return
        total = sum(dt for _, dt in self.records)
        logger.info("profile summary: total %.3fs across %d stage(s)", total, len(self.records))
        for label, dt in self.records:
            logger.info("profile: %-34s %.3fs", label, dt)


def _meta_identity_fingerprint(meta):
    """
    Cheap detect of Table.meta mutations (top-level key identity).

    Omits ``mask``: it is recomputed after ``load_metadata`` each run, so its
    object id would otherwise force a cache write every iteration.
    """
    skip = {"mask"}
    return tuple(sorted((str(k), id(v)) for k, v in meta.items() if k not in skip))


def do_trim_metadata(sndata, keep_keys=TRIM_METADATA_KEEP_KEYS):
    """
    Remove metadata entries whose object name is not in the table.

    Parameters
    ----------
    sndata : dict of Table
        Per-sheet tables with .meta dicts.
    keep_keys : set or list, optional
        Metadata keys to always keep (e.g. 'all_sndata').

    Returns
    -------
    dict
        sndata (modified in place).
    """
    for key in sndata.keys():
        init_keys = list(sndata[key].meta.keys())
        removed = []
        for obj in init_keys:
            if obj in keep_keys:
                continue
            if obj not in sndata[key]['Name'].data:
                removed.append(obj)
                del sndata[key].meta[obj]
        if removed:
            tail = removed[:8]
            extra = f" (+{len(removed) - len(tail)} more)" if len(removed) > len(tail) else ""
            logger.info("trim %r: removed %d stale meta key(s): %s%s",
                        key, len(removed), ", ".join(tail), extra)
    return sndata


def main(
    redo=False,
    download_yse=True,
    update_classification=False,
    alert=True,
    yse_sql_query=DEFAULT_YSE_SQL_QUERY,
    always_update=False,
    trim_metadata=False,
    redo_obj=None,
    metadata=DEFAULT_METADATA,
    profile=False,
):
    """
    Run the full progenitor pipeline: Sheets, metadata, optional classification and alerts.

    Parameters
    ----------
    redo : bool, optional
        Re-fetch metadata even if already present.
    download_yse : bool, optional
        Add YSE targets to sheet data.
    update_classification : bool, optional
        Reclassify objects across sheets and upload. When ``True``, the TNS
        public-objects CSV zip is downloaded once at the start; after the main
        metadata pass, each transient is crossmatched to that catalog by name only
        (TNS ``name`` / ``internal_names``) and the sheet ``Classification``
        column is set from the
        TNS ``type`` field before ``gather_type`` runs.
    alert : bool, optional
        When ``update_classification`` is True, send email for progenitor-interest
        transients that **newly appear** on one of the core science tabs (Type Ia,
        Ib/c, II-P/II-L, IIn, IIb) compared to ``meta['curr_table']``.
    yse_sql_query : str, optional
        SQL for YSE target query.
    always_update : bool, optional
        Always save metadata after each type.
    trim_metadata : bool, optional
        Call do_trim_metadata before final save.
    redo_obj : list of str, optional
        Object names to redo metadata for only.
    metadata : list of str, optional
        Metadata types to fetch (e.g. 'ned', 'distance', 'tns').
    profile : bool, optional
        Log per-stage wall-clock timings.
    """
    profiler = _Profiler(enabled=profile)
    tns_public_table = None
    if update_classification:
        profiler.start("tns_catalog_download")
        logger.info("Downloading TNS public objects catalog for reclassification…")
        tns_public_table = util.fetch_tns_public_objects_table()
        logger.info(
            "TNS public catalog loaded: %d object(s)", len(tns_public_table)
        )
        profiler.stop("tns_catalog_download", "tns catalog download")

    profiler.start("download_progenitor_data")
    sndata = sheetproc.download_progenitor_data(util.params['SHEET'])
    profiler.stop("download_progenitor_data", "download progenitor data")
    profiler.start("load_metadata")
    sndata = sheetproc.load_metadata_from_file(sndata, util.params['metadata'])
    profiler.stop("load_metadata", "load metadata from cache")

    for key in sndata.keys():
        if 'all_sndata' in sndata[key].meta.keys():
            meta = copy.copy(sndata[key].meta)
            curr_table = Table(sndata[key], meta={})
            all_sn_data = Table(sndata[key].meta['all_sndata'], meta={})
            sndata[key] = all_sn_data
            sndata[key].meta = meta
            sndata[key].meta['curr_table'] = curr_table

    if download_yse:
        profiler.start("add_yse_targets")
        sndata = util.add_yse_targets(sndata, yse_sql_query=yse_sql_query)
        profiler.stop("add_yse_targets", "add yse targets")

    all_keys = list(sndata.keys())
    defer_upload_until_reclassification = update_classification
    _meta_count_skip = set(TRIM_METADATA_KEEP_KEYS) | {"name"}
    for type_key in all_keys:
        n_meta_buckets = sum(
            1 for k in sndata[type_key].meta if k not in _meta_count_skip
        )
        logger.info(
            "sheet %r: ~%d object meta bucket(s); metadata passes: %s",
            type_key,
            n_meta_buckets,
            ",".join(metadata),
        )

        for dattype in metadata:
            step_key = f"{type_key}:metadata:{dattype}"
            profiler.start(step_key)
            update = util.add_metadata(sndata[type_key], dattype, redo=redo, redo_obj=redo_obj)
            profiler.stop(step_key, f"{type_key}: metadata {dattype}")
            if vars(update) != vars(sndata[type_key]) or always_update:
                sndata[type_key] = copy.copy(update)
                save_key = f"{type_key}:save_meta_loop"
                profiler.start(save_key)
                sheetproc.save_metadata_to_file({type_key: sndata[type_key]}, util.params['metadata'])
                profiler.stop(save_key, f"{type_key}: save metadata cache (loop)")

        cols = sheetproc.params['cols']
        logger.info("sheet %r: column pipeline (%d): %s", type_key, len(cols), ",".join(cols))
        fp_before_cols = _meta_identity_fingerprint(sndata[type_key].meta)
        profiler.start(f"{type_key}:column_pipeline")
        for coltype in cols:
            sndata[type_key] = util.add_data(sndata[type_key], coltype)
        profiler.stop(f"{type_key}:column_pipeline", f"{type_key}: column pipeline")

        if len(sndata[type_key]) > 1:
            sndata[type_key] = unique(sndata[type_key], keys='Name')
        sndata[type_key].sort('Discovery Date')
        fp_after_cols = _meta_identity_fingerprint(sndata[type_key].meta)
        sndata[type_key].meta['mask'] = util.get_classification_mask(sndata[type_key])

        # For reclassification runs, this intermediate upload is overwritten by
        # the final all-tab upload below; defer it to avoid duplicate Sheets I/O.
        if not defer_upload_until_reclassification:
            profiler.start(f"{type_key}:upload_intermediate")
            # Upload only this tab — full-workbook upload each iteration exceeds the
            # default 60 Sheets write requests/minute/user quota.
            sheetproc.upload_progenitor_data(
                util.params['SHEET'], {type_key: sndata[type_key]}, mask=True
            )
            profiler.stop(f"{type_key}:upload_intermediate", f"{type_key}: upload intermediate tab")
        if fp_after_cols != fp_before_cols:
            logger.info(
                "sheet %r: meta changed after column pipeline; saving cache",
                type_key,
            )
            profiler.start(f"{type_key}:save_meta_postcols")
            sheetproc.save_metadata_to_file({type_key: sndata[type_key]}, util.params['metadata'])
            profiler.stop(f"{type_key}:save_meta_postcols", f"{type_key}: save metadata cache (post-cols)")

    if trim_metadata:
        profiler.start("trim_metadata")
        do_trim_metadata(sndata)
        profiler.stop("trim_metadata", "trim metadata")
    profiler.start("save_metadata_all")
    sheetproc.save_metadata_to_file(sndata, util.params['metadata'], make_backup=True)
    profiler.stop("save_metadata_all", "save metadata all (+backup)")

    if update_classification:
        if tns_public_table is not None and len(tns_public_table) > 0:
            profiler.start("apply_tns_catalog")
            util.apply_tns_public_catalog_to_sndata(sndata, tns_public_table)
            profiler.stop("apply_tns_catalog", "apply tns public catalog")
        else:
            logger.warning(
                "Skipping TNS public-catalog crossmatch (no catalog loaded)"
            )

        new_sndata = {}
        for key in all_keys:
            profiler.start(f"reclassify:{key}")
            new_sndata[key] = util.gather_type(sndata, key)
            profiler.stop(f"reclassify:{key}", f"reclassify {key}")
            logger.info("reclassify %r: %d row(s)", key, len(new_sndata[key]))

        already_placed = set()
        for k in all_keys:
            already_placed.update(str(x) for x in new_sndata[k]["Name"].data)
        for key in all_keys:
            for row in sndata[key]:
                nm = str(row["Name"])
                if nm in already_placed:
                    continue
                new_sndata["Other"].add_row(row)
                new_sndata["Other"].meta[row["Name"]] = sndata[key].meta[row["Name"]]
                already_placed.add(nm)

        for key in new_sndata.keys():
            new_sndata[key] = unique(new_sndata[key], keys='Name')
            new_sndata[key].sort('Discovery Date')
            new_sndata[key].meta['mask'] = util.get_classification_mask(new_sndata[key])
            if (
                alert
                and key in ALERT_RECLASSIFICATION_SHEET_KEYS
                and "curr_table" in new_sndata[key].meta.keys()
            ):
                mask = new_sndata[key].meta["mask"]
                curr_table = new_sndata[key].meta["curr_table"]
                curr_names_set = set(curr_table["Name"].data)
                for ii, row in enumerate(new_sndata[key]):
                    if not mask[ii] and row["Name"] not in curr_names_set:
                        util.post_alert(row)

        profiler.start("upload_reclassified")
        sheetproc.upload_progenitor_data(util.params['SHEET'], new_sndata, mask=True)
        profiler.stop("upload_reclassified", "upload reclassified sheets")
        for key in new_sndata.keys():
            copy_table_data = copy.copy(new_sndata[key])
            copy_table_data.meta = None
            new_sndata[key].meta['all_sndata'] = copy_table_data
        profiler.start("save_metadata_reclassified")
        sheetproc.save_metadata_to_file(new_sndata, util.params['metadata'], make_backup=True)
        profiler.stop("save_metadata_reclassified", "save metadata reclassified (+backup)")
    profiler.log_summary()
