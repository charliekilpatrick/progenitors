import argparse
import sys
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import units as u
from astropy.time import Time

def message(msg):
    print('\n\n'+msg+'\n'+'#'*80+'\n'+'#'*80+'\n\n')

def parse_two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    values = map(float, values)
    return(values)

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in ra and ':' in dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('--redo', default=False, action='store_true',
                        help='Redo metadata for all objects in database.')
    parser.add_argument('--alert', default=False, action='store_true',
                        help='Send an email alert if there are new YSE objs.')
    parser.add_argument('--always-update', default=False, action='store_true',
                        help='Send an email alert if there are new YSE objs.')
    parser.add_argument('--update-classification', default=False,
                        action='store_true',
                        help='Update the classification of objects.')
    parser.add_argument('--yse-sql-query', default='378', type=str,
                        help='YSE SQL query to use for grabbing YSE data.')
    parser.add_argument('--trim-metadata', default=False, action='store_true',
                        help='Trim keys from metadata that are not in table.')
    parser.add_argument('--redo-obj', default=None, type=str,
                        help='Redo metadata for these objects.  Input should '+\
                        'be a comma-separated list of objs (1987A,1993J).')
    parser.add_argument('--metadata', default=None, type=str,
                        help='Only do metadata for these data types.')

    args = parser.parse_args()

    return(args)
