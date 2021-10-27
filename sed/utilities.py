from astropy.coordinates import SkyCoord
from astropy.table import Column
from astropy import units as u
import numpy as np
import dustmaps.sfd
import os


def import_dustmap():
    for pole in ['ngp', 'sgp']:
        datadir=dustmaps.sfd.data_dir()
        file = os.path.join(datadir, 'sfd',
            'SFD_dust_4096_{}.fits'.format(pole))
        if not os.path.exists(file):
            dustmaps.sfd.fetch()
            break

def get_sfd():
    return(dustmaps.sfd.SFDQuery())

# Round a floating point number f to n significant digits
def round_to_n(f, n):
    code = '%.{0}g'.format(int(n))
    val = float(code % f)
    return(val)

# For an arbitrary input Ra/Dec, output a SkyCoord object
def parse_coord(ra, dec):
    def check_coord(ra, dec):
        if ':' in str(ra) and ':' in str(dec):
            coord = SkyCoord(ra, dec, unit=(u.hour, u.deg))
            return(coord)
        else:
            try:
                ra = float(ra) ; dec = float(dec)
                coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
                return(coord)
            except:
                return(None)

    if (isinstance(ra, (list, np.ndarray, Column)) and
        isinstance(dec, (list, np.ndarray, Column))):
        coords = np.array([check_coord(r,d) for r, d in zip(ra, dec)])
        return(coords)
    else:
        return(check_coord(ra, dec))
