"""
Generate a generic HR diagram from the package progenitors.dat (no input file required).

Run from repository root:
    python -m progenitors.sed.plotting.plot_hr_progenitors

Output: progenitors/sed/figures/progenitors_hr.png
"""
from .hr import plot_hr_from_progenitors


def main():
    outpath = plot_hr_from_progenitors()
    print("Saved:", outpath)


if __name__ == "__main__":
    main()
