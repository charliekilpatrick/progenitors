"""Tests for pipeline profiling mode."""
import copy

from astropy.table import Table


def _mk_table(name: str) -> Table:
    t = Table(rows=[(name, "00:00", "+00", "2020-01-01")], names=("Name", "RA", "Dec", "Discovery Date"))
    t.meta = {name: {"ok": True}, "name": name}
    return t


def test_pipeline_profile_logs_summary(monkeypatch, caplog):
    from progenitors import pipeline

    data = {"Type Ia": _mk_table("SN-Ia-1"), "Other": _mk_table("SN-Other-1")}

    monkeypatch.setattr(pipeline.sheetproc, "download_progenitor_data", lambda _sid: copy.deepcopy(data))
    monkeypatch.setattr(pipeline.sheetproc, "load_metadata_from_file", lambda d, _dir: d)
    monkeypatch.setattr(pipeline.sheetproc, "save_metadata_to_file", lambda *_a, **_k: None)
    monkeypatch.setattr(pipeline.util, "add_yse_targets", lambda d, yse_sql_query=None: d)
    monkeypatch.setattr(pipeline.util, "add_metadata", lambda d, *_a, **_k: d)
    monkeypatch.setattr(pipeline.util, "add_data", lambda d, *_a, **_k: d)
    monkeypatch.setattr(pipeline.util, "get_classification_mask", lambda d: [False] * len(d))
    monkeypatch.setitem(pipeline.sheetproc.params, "cols", [])
    monkeypatch.setattr(pipeline.sheetproc, "upload_progenitor_data", lambda *_a, **_k: None)

    caplog.set_level("INFO")
    pipeline.main(
        download_yse=False,
        update_classification=False,
        metadata=[],
        profile=True,
    )

    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "profile summary:" in text
    assert "profile: download progenitor data" in text
    assert "profile: load metadata from cache" in text
