import requests
import numpy as np
import dustmaps.sfd
import pandas
import json
import os
import re
from collections import OrderedDict
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy.time import Time
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.io import ascii

from astroquery.ned import Ned
from astroquery.vizier import Vizier

Vizier.ROW_LIMIT = -1

import warnings
warnings.filterwarnings('ignore')

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

def get_coord(row):

    coord = None
    if 'RA' in row.keys() and 'Dec' in row.keys():
        coord = parse_coord(row['RA'], row['Dec'])

    return(coord)

def check_glade(row, table=None):
    coord = get_coord(row)
    if not coord:
        return({'error': 'no coordinate'})

    result = Vizier.query_region(coord, radius=60*u.arcmin,
        catalog='VII/281/glade2')

    if len(result)==0:
        return({'error': 'no result'})

    result = result[0]
    # Get rid of columns that we won't use
    result = result['GWGC','HyperLEDA','z','RAJ2000','DEJ2000','Dist','BMAG']

    # Add objects because GLADE does not contain M31, M33
    result.add_row(['M31','M31',-0.001001,10.6847929,41.2690650,0.731,-20.10])
    result.add_row(['M33','M33',-0.000597,23.4620417,30.6602222,0.832,-18.31])

    mask = ~result['GWGC'].mask & ~result['z'].mask & (result['GWGC']!='---')
    result = result[mask]

    sep = []
    # Approximate separation in kpc
    distmask = result['Dist'].mask
    for i,row in enumerate(result):
        c1 = SkyCoord(row['RAJ2000'],row['DEJ2000'],unit='deg')
        if not distmask[i]:
            dist = row['Dist'] * 1000.0
        # Ignore rows that have negative redshift - probably not host anyway
        elif row['z']<0.0:
            dist = 9.0e9
        else:
            dist = row['z']*100.0*43.4*1000.0
        sep.append(c1.separation(coord).radian * dist)

    result.add_column(Column(sep, name='Separation'))

    # Impose a cut at a projected separation of 50 kpc
    result = result[result['Separation']<50.0]

    if len(result)==0:
        return({'error': 'no candidates'})

    # Try easy case where there's only one obvious candidate
    if len(result)==1:
        best = result[0]
        url='https://ned.ipac.caltech.edu/?q=byname&objname={0}'
        url=url.format(best['GWGC'].replace(' ','+'))
        output = {'name': best['GWGC'], 'ra': best['RAJ2000'],
                  'DEC': best['DEJ2000'], 'redshift': best['z'],
                  'error': None, 'distance': best['Dist'],
                  'separation': best['Separation'], 'url': url}
        return(output)

    mask = ~result['GWGC'].mask & ~result['z'].mask & ~result['BMAG'].mask &\
        (result['GWGC']!='---')
    result = result[mask]
    if len(result)==0:
        return({'error': 'no candidates'})

    Lb = 3.839e33

    # For optimization statistic, use B-band luminosity / separation**2
    newcol = Column(Lb*10**(-0.4*(result['BMAG']-5.48))/result['Separation']**2,
        name='Surface Brightness')
    result.add_column(newcol)

    best = np.argmax(result['Surface Brightness'])
    best = result[best]

    url='https://ned.ipac.caltech.edu/?q=byname&objname={0}'
    url=url.format(best['GWGC'].replace(' ','+'))
    output = {'name': best['GWGC'], 'ra': best['RAJ2000'],
              'DEC': best['DEJ2000'], 'redshift': best['z'],
              'error': None, 'distance': best['Dist'],
              'separation': best['Separation'], 'url': url}

    return(output)


row = {'RA': '10:16:56.667', 'Dec': '+73:23:51.290'}

output = check_glade(row)
print(output)
