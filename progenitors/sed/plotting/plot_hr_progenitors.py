"""
Generate a generic HR diagram from the package progenitors.dat (no input file required).

Run from repository root:
    python -m progenitors.sed.plotting.plot_hr_progenitors

MIST single-star tracks are read from ``progenitors/sed/data/mist/`` (see ``hr.MIST_MASSES``).
Set ``MIST_DIR`` only if your grid lives elsewhere. Output: ``progenitors/sed/figures/progenitors_hr.eps``.
"""
import os

from .hr import plot_hr_from_progenitors, _sed_data_dir


def main():
    mist_dir = os.path.normpath(os.path.join(_sed_data_dir(), "mist"))
    outpath = plot_hr_from_progenitors(mist_dir=mist_dir)
    print("Saved:", outpath)


if __name__ == "__main__":
    main()
