"""Performance-oriented behavior tests for pipeline upload flow."""
import copy

from astropy.table import Table


def _mk_table(name: str) -> Table:
    t = Table(rows=[(name, "00:00", "+00", "2020-01-01")], names=("Name", "RA", "Dec", "Discovery Date"))
    t.meta = {name: {"ok": True}, "name": name}
    return t


def test_update_classification_defers_intermediate_uploads(monkeypatch):
    """Classification run should upload once (final), not once per tab + final."""
    from progenitors import pipeline

    data = {
        "Type Ia": _mk_table("SN-Ia-1"),
        "Other": _mk_table("SN-Other-1"),
    }

    monkeypatch.setattr(pipeline.util, "fetch_tns_public_objects_table", lambda: Table())
    monkeypatch.setattr(pipeline.sheetproc, "download_progenitor_data", lambda _sid: copy.deepcopy(data))
    monkeypatch.setattr(pipeline.sheetproc, "load_metadata_from_file", lambda d, _dir: d)
    monkeypatch.setattr(pipeline.sheetproc, "save_metadata_to_file", lambda *_a, **_k: None)
    monkeypatch.setattr(pipeline.util, "add_yse_targets", lambda d, yse_sql_query=None: d)
    monkeypatch.setattr(pipeline.util, "add_metadata", lambda d, *_a, **_k: d)
    monkeypatch.setattr(pipeline.util, "add_data", lambda d, *_a, **_k: d)
    monkeypatch.setattr(pipeline.util, "get_classification_mask", lambda d: [False] * len(d))
    monkeypatch.setattr(pipeline.util, "gather_type", lambda snd, key: Table(snd[key], meta=copy.copy(snd[key].meta)))
    monkeypatch.setattr(pipeline.util, "apply_tns_public_catalog_to_sndata", lambda *_a, **_k: None)
    monkeypatch.setitem(pipeline.sheetproc.params, "cols", [])

    calls = []

    def _upload(_sid, all_tables, mask=False):
        calls.append((set(all_tables.keys()), mask))

    monkeypatch.setattr(pipeline.sheetproc, "upload_progenitor_data", _upload)

    pipeline.main(
        redo=False,
        download_yse=False,
        update_classification=True,
        alert=False,
        always_update=False,
        trim_metadata=False,
        metadata=[],
    )

    assert len(calls) == 1
    keys, mask = calls[0]
    assert keys == set(data.keys())
    assert mask is True
