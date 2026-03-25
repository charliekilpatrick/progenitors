"""SFD dust map availability (shared by pipeline util and SED code)."""
import os


def import_dustmap():
    """
    Ensure SFD dust map data is available; fetch if missing.

    Checks for SFD FITS files in the dustmaps data directory and calls
    ``dustmaps.sfd.fetch()`` if needed.
    """
    import dustmaps.sfd

    for pole in ['ngp', 'sgp']:
        datadir = dustmaps.sfd.data_dir()
        path = os.path.join(datadir, 'sfd',
                            'SFD_dust_4096_{}.fits'.format(pole))
        if not os.path.exists(path):
            dustmaps.sfd.fetch()
            break


def get_sfd():
    """Return an ``SFDQuery`` instance for E(B-V) lookups."""
    import dustmaps.sfd

    return dustmaps.sfd.SFDQuery()
