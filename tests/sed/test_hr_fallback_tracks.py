"""MIST HR fallback when CMD grid is not installed."""
from progenitors.sed.plotting import hr


def test_load_mist_tracks_fallback_without_mist_directory(tmp_path, monkeypatch):
    monkeypatch.setenv("MIST_DIR", str(tmp_path / "no_mist_here"))
    t = hr._load_mist_tracks()
    assert t is not None
    assert len(t) > 0
    for col in ("star_age", "log_Teff", "log_L", "mass"):
        assert col in t.colnames
    masses = sorted(set(float(m) for m in t["mass"]))
    assert masses == [float(m) for m in hr.MIST_MASSES]
