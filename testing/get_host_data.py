#import util
from astropy.table import Table
import requests
import numpy as np
import pandas

from astroquery.ned import Ned

row = {'Name': 'test'}

# Round a floating point number f to n significant digits
def round_to_n(f, n):
    code = '%7.{0}f'.format(int(n))
    val = code % float(f)
    return(val.strip())

dm_hierarchy = [{'name':'Cepheids','keys':['Cepheids','Cepheid']},
                {'name':'TRGB','keys':['TRGB']},
                {'name':'Tully-Fisher','keys':['Tully-Fisher','Tully est']},
                {'name':'Fundamental Plane','keys':['FP']},
                {'name':'SNIa','keys':['SNIa']},
                {'name':'SNII','keys':['SNII optical','SNII radio']},
                {'name':'Sosies','keys':['Sosies']},
                {'name':'Planetary Nebula Luminosity Function',
                    'keys':['PNLF']},
                {'name':'Ring Diameter','keys':['Ring Diameter']}]

def get_host_distance(name):

    df = []
    if name is None or not name:
        name = row['Host']

    # If there is not host data then can't do a search for NED data
    if name is None or not name or name in ['---']:
        return({})

    url = 'https://ned.ipac.caltech.edu/cgi-bin/nDistance?name={0}'
    name_fmt = name.strip().replace(' ','+')
    url = url.format(name_fmt)

    try:
        r = requests.get(url)
    except requests.exceptions.ConnectionError:
        print('Connection error for: {0}'.format(row['Name']))
        return({})

    if r.status_code==200:
        try:
            df = pandas.read_html(r.content)
        except ValueError:
            pass

    if len(df)>1:
        table = Table.from_pandas(df[1])
        return({'table': table})
    else:
        return({'table': None})

    return({})

def get_best_distance(table):

    if table:
        mask = table['err(m-M)']>0
        table = table[mask]

        for m in dm_hierarchy:
            if any([t in table['Method'] for t in m['keys']]):
                mask = np.array([l in m['keys'] for l in table['Method']])
                best = table[mask]
                idx = np.argmax([row['REFCODE'][0:4] for row in best])
                best = best[idx]

                dist = 10**(0.2*(best['(m-M)']-25.0))
                dmerr = best['err(m-M)']
                disterr = 4.6052e-6*np.exp(0.46052*best['(m-M)'])*dmerr

                output = {'distance': round_to_n(dist, 2),
                    'distance_error': round_to_n(disterr, 2),
                    'method': m['name'], 'reference': best['REFCODE']}

                return(output)

    if (row['Redshift'] and is_number(row['Redshift']) and
        float(row['Redshift'])>0):

        dist = cosmo.luminosity_distance(float(row['Redshift'])).value

        # error is dominated by H0 with H0=67.8+/-0.9
        # https://arxiv.org/pdf/1502.01589.pdf
        unc_factor = 0.9/67.8
        disterr = dist * unc_factor

        # Add systematic uncertainty from Scolnic+2018
        # https://arxiv.org/pdf/1710.00845.pdf
        sys_velo = 250.0/2.998e5 # km/s
        sys_unc_factor = sys_velo / float(row['Redshift'])
        sys_disterr = dist * sys_unc_factor

        disterr = np.sqrt(disterr**2 + sys_disterr**2)

        output = {'distance': round_to_n(dist, 2),
                'distance_error': round_to_n(disterr, 2),
                'method': 'Redshift', 'reference': '2016A%26A...594A..13P'}

        return(output)

    return({})

for obj in ['NGC4594']:

    data = Ned.query_object(obj)
    ra = data['RA'][0]
    dec = data['DEC'][0]

    table = get_host_distance(obj)
    dist = get_best_distance(table['table'])

    distance = dist['distance']
    dist_err = dist['distance_error']

    url='https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname={0}'.format(obj)
    print(obj, ra, dec, url, f'{distance} ({dist_err})', dist['method'], dist['reference'])
