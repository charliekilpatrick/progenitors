"""Unit tests for sheetproc module."""
import os
import pytest

# Allow tests to run without real credentials
if "PROGENITORS_SHEET" not in os.environ:
    os.environ.setdefault("PROGENITORS_SHEET", "test-dummy")
if "YSE_USER" not in os.environ:
    os.environ.setdefault("YSE_USER", "test")
if "YSE_PASSWORD" not in os.environ:
    os.environ.setdefault("YSE_PASSWORD", "test")


def test_transpose():
    from progenitors.sheets.sheetproc import transpose
    data = [[1, 2, 3], [4, 5, 6]]
    out = transpose(data)
    assert out == [[1, 4], [2, 5], [3, 6]]


def test_params_cols():
    from progenitors.sheets.sheetproc import params
    assert "cols" in params
    assert "Name" in params["cols"]
    assert "SHEET" in params


def test_convert_table_to_lists():
    from astropy.table import Table
    from progenitors.sheets.sheetproc import convert_table_to_lists, params
    t = Table({"Name": ["SN1", "SN2"], "RA": ["01:00:00", "02:00:00"], "Dec": ["+00:00", "+01:00"]})
    out = convert_table_to_lists(t)
    assert out[0] == list(params["cols"])  # header row
    assert len(out) == 3  # header + 2 rows
    assert "SN1" in str(out[1])


def test_load_save_metadata_roundtrip(tmp_path):
    from astropy.table import Table
    from progenitors.sheets.sheetproc import (
        METADATA_FILE_SUFFIX,
        load_metadata_from_file,
        save_metadata_to_file,
    )

    all_data = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    all_data["Type Ia"].meta = {"SN1": {"test": 1}}
    save_metadata_to_file(all_data, str(tmp_path))
    zst = tmp_path / f"Type_Ia{METADATA_FILE_SUFFIX}"
    assert zst.is_file(), f"expected {zst}"
    loaded = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    loaded = load_metadata_from_file(loaded, str(tmp_path))
    assert loaded["Type Ia"].meta.get("SN1", {}).get("test") == 1


def test_load_metadata_recover_from_backup_when_main_corrupt(tmp_path, monkeypatch):
    """Corrupt primary .pkl.zst should fall back to backup/*.backup.pkl.zst."""
    from astropy.table import Table

    from progenitors.sheets.sheetproc import (
        METADATA_BACKUP_SUFFIX,
        METADATA_FILE_SUFFIX,
        load_metadata_from_file,
        _atomic_metadata_zstd_dump,
    )

    monkeypatch.setenv("PROGENITORS_METADATA_QUARANTINE_CORRUPT", "1")
    savekey = "Type_Ia"
    backup_dir = tmp_path / "backup"
    backup_dir.mkdir()
    main = tmp_path / f"{savekey}{METADATA_FILE_SUFFIX}"
    backup = backup_dir / f"{savekey}{METADATA_BACKUP_SUFFIX}"
    meta_good = {"SN1": {"recovered": True}}
    _atomic_metadata_zstd_dump(meta_good, str(backup))
    main.write_bytes(b"\x00\x01\x02not-zstd")

    all_data = {"Type Ia": Table({"Name": ["SN1"], "RA": ["00:00"], "Dec": ["+00"]})}
    loaded = load_metadata_from_file(all_data, str(tmp_path))
    assert loaded["Type Ia"].meta["SN1"]["recovered"] is True
    assert (tmp_path / f"{savekey}{METADATA_FILE_SUFFIX}.corrupt").is_file()


def test_load_metadata_legacy_pkl(tmp_path):
    """Uncompressed .pkl is still readable until migrated by a save."""
    import pickle

    from astropy.table import Table
    from progenitors.sheets.sheetproc import (
        LEGACY_METADATA_SUFFIX,
        METADATA_FILE_SUFFIX,
        load_metadata_from_file,
        save_metadata_to_file,
    )

    pkl = tmp_path / f"Type_Ia{LEGACY_METADATA_SUFFIX}"
    with open(pkl, "wb") as f:
        pickle.dump({"SN9": {"legacy": True}}, f)

    loaded = {"Type Ia": Table({"Name": ["SN9"], "RA": ["00:00"], "Dec": ["+00"]})}
    loaded = load_metadata_from_file(loaded, str(tmp_path))
    assert loaded["Type Ia"].meta["SN9"]["legacy"] is True

    save_metadata_to_file(loaded, str(tmp_path))
    assert (tmp_path / f"Type_Ia{METADATA_FILE_SUFFIX}").is_file()
    assert not pkl.is_file()


def test_load_save_metadata_parallel_sheets(tmp_path):
    """Multiple sheet keys round-trip (exercises threaded load/save paths)."""
    from astropy.table import Table
    from progenitors.sheets.sheetproc import load_metadata_from_file, save_metadata_to_file

    all_data = {
        "Type Ia": Table({"Name": ["A"], "RA": ["00:00"], "Dec": ["+00"]}),
        "Type II": Table({"Name": ["B"], "RA": ["01:00"], "Dec": ["+01"]}),
    }
    all_data["Type Ia"].meta = {"A": {"n": 1}}
    all_data["Type II"].meta = {"B": {"n": 2}}
    save_metadata_to_file(all_data, str(tmp_path))

    loaded = {
        "Type Ia": Table({"Name": ["A"], "RA": ["00:00"], "Dec": ["+00"]}),
        "Type II": Table({"Name": ["B"], "RA": ["01:00"], "Dec": ["+01"]}),
    }
    loaded = load_metadata_from_file(loaded, str(tmp_path))
    assert loaded["Type Ia"].meta["A"]["n"] == 1
    assert loaded["Type II"].meta["B"]["n"] == 2


def test_atomic_pickle_dump_keeps_existing_on_failure(tmp_path, monkeypatch):
    """Failed serialize must not truncate an existing target file."""
    import pickle

    from progenitors.sheets import sheetproc as sp

    path = tmp_path / "existing.pkl"
    path.write_bytes(b"KEEP_ME")

    def _fail_dump(*_args, **_kwargs):
        raise RuntimeError("simulated write failure")

    monkeypatch.setattr(pickle, "dump", _fail_dump)
    with pytest.raises(RuntimeError):
        sp._atomic_pickle_dump({"x": 1}, str(path))
    assert path.read_bytes() == b"KEEP_ME"
