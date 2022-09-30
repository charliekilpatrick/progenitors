import requests
import numpy as np
import dustmaps.sfd
import pandas
import json
import os
import re
from collections import OrderedDict
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, unique
from astropy.time import Time
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.io import ascii

import astroquery
from astroquery.mast import Observations
from astroquery.ned import Ned
from astroquery.vizier import Vizier

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import smtplib

import warnings
warnings.filterwarnings('ignore')

Ned.TIMEOUT=30
Vizier.ROW_LIMIT = -1

basedir='/home/ckilpatrick/scripts/python/progenitors/'
params = {
    'SHEET':'1paDfeYsJyv9X_XL26gV9Pk2xF9VjG70T4Wly5lRPuxs',
    'SHEET_TEST': '1ZTASxDCwj7H_rcjK1rU4Vi-FnnNkiaU5hoMggr0pn3o',
    'token':'/home/ckilpatrick/scripts/python/progenitors/token.pickle',
    'target': '/home/ckilpatrick/scripts/python/progenitors/',
    'tns': {'api_key': 'eaa6718c10187b6980f52bac3d79e8f5e9ba8e0f'},
    'yse': {'user': 'ckilpatrick', 'password': 'Vfg190OW@K9E*g4$Bpmw'},
    'metadata': basedir+'metadata/',
    'extinction': basedir+'extinction/',
}

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

authcode='GGLWYREIaQMBDac8dfrKKm5KIdm3yj88IEGclIaL'
ads_headers={'Authorization': 'Bearer '+authcode}

email_args = {
    'from_addr': 'Supernova Progenitor Alerts',
    'smtpserver': '%s:%s' % ('smtp.gmail.com', 587),
    'gmail_login': 'hst.supernovae@gmail.com',
    'gmail_password': 'okwqoshxzbokebxz',
    'to_addr': 'ckilpatrick@northwestern.edu',
    'subject': 'Supernova Progenitor Target Summary'
}

email_msg = '''<html><body><p>Bright transients for APF follow up</p>
{targets}
<p>CDK</p>
</body></html>'''

def sendEmail(from_addr, to_addr, subj, message, login, password, smtpserver):

    print('Preparing email')

    msg = MIMEMultipart('alternative')
    msg['Subject'] = subj
    msg['From'] = from_addr
    msg['To'] = to_addr
    payload = MIMEText(message, 'html')
    msg.attach(payload)

    with smtplib.SMTP(smtpserver) as server:
        try:
            server.starttls()
            server.login(login, password)
            resp = server.sendmail(from_addr, [to_addr], msg.as_string())
            print('Send email success: {0}'.format(to_addr))
        except Exception as e:
            print('Send email fail: {0}'.format(to_addr))
            print('ERROR:',e)

    return(1)

def is_number(val):
    try:
        val = float(val)
        return(True)
    except (ValueError, TypeError):
        return(False)

def check_dict(dic, keys):

    for key in keys:
        if isinstance(dic, dict) and key in dic.keys():
            dic = dic[key]
        else:
            return(None)

    return(dic)

def search_tns(coord, radius=5):
    ra, dec = coord.to_string(style='hmsdms', sep=':', precision=3).split()
    data = {'ra': ra, 'dec': dec, 'radius': str(radius)}

    json_file=OrderedDict(data)
    post_data=[('api_key',(None, params['tns']['api_key'])),
                      ('data',(None,json.dumps(json_file)))]

    url='https://www.wis-tns.org/api/get/search'

    response=requests.post(url, files=post_data)

    if response.status_code==200:
        out=json.loads(response.content)
        if len(out['data']['reply'])>0:
            return(out['data']['reply'][0]['objname'])

    return('')

def get_tns_data(row, table=None):
    USER_ID='97993'
    USER_NAME='YSE_progenitors'
    data = {'objname': row['Name'], 'spectra': '1', 'photometry': '1'}
    headers={'User-Agent':'tns_marker{"tns_id":'+str(USER_ID)+\
        ', "type":"bot", "name":"'+USER_NAME+'"}'}
    json_file=OrderedDict(data)
    post_data=[('api_key',(None, params['tns']['api_key'])),
               ('data',(None,json.dumps(json_file)))]
    url='https://www.wis-tns.org/api/get/object'
    r=requests.post(url, files=post_data, headers=headers)
    if r.status_code==200:
        out = json.loads(r.content)
        out = check_dict(out, ['data','reply'])

        if not out:
            return({'error': 'bad dict'})

        disc = get_tns_source(row, 'discovery')
        clas = get_tns_source(row, 'classification')

        if disc: out['discovery_reference']=disc
        if clas: out['classification_reference']=clas

        objpage = 'https://www.wis-tns.org/object/{0}'
        out['url']=objpage.format(row['Name'])

        return(out)
    else:
        print(r.status_code)
        print(r.content)

    return({'error': 'timeout data'})

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
    code = '%7.{0}f'.format(int(n))
    val = code % float(f)
    return(val.strip())

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
    if 'RA' in row.colnames and 'Dec' in row.colnames:
        coord = parse_coord(row['RA'], row['Dec'])

    return(coord)

def get_osc(row, table=None):

    url = check_dict(table.meta, [row['Name'],'osc','url'])
    if url: return(url)

    return('')

def get_tns(row, table=None):

    url = check_dict(table.meta, [row['Name'],'tns','url'])
    if url: return(url)

    return('')

def get_yse(row, table=None):

    url = check_dict(table.meta, [row['Name'],'yse','url'])
    if url: return(url)

    return('')

def get_yse_data(row, table=None):
    name = row['Name']
    user = params['yse']['user']
    password = params['yse']['password']

    url = 'https://ziggy.ucolick.org/yse/download_data/{0}'
    url = url.format(name.strip())

    page = url.replace('download_data', 'transient_detail')

    out = {}
    try:
        r = requests.get(page, auth=(user, password), timeout=180)
        if r.status_code==200:
            out['url']=page
    except requests.exceptions.ReadTimeout:
        print('Timeout on YSE page for: {0}'.format(row['Name']))
        return({'error': 'timeout data'})
    except requests.exceptions.ConnectionError:
        print('Connection error on YSE page for: {0}'.format(row['Name']))
        return({'error': 'timeout data'})

    return(out)

    try:
        r = requests.get(url, auth=(user, password), timeout=90)
    except requests.exceptions.ReadTimeout:
        out['error']='timeout data'
        print('Timeout on YSE data for: {0}'.format(row['Name']))
        return(out)
    except requests.exceptions.ConnectionError:
        out['error']='timeout data'
        print('Connection error on YSE page for: {0}'.format(row['Name']))
        return(out)

    if r.status_code==200:
        data = json.loads(r.content)
        out.update(data[name])
        out['download']=url
        return(out)

    return({'error': r.status_code})

def add_data(table, add_key):

    if add_key not in table.keys():
        newcol = []
        for row in table:
            newcol.append(get_data[add_key](row, table=table))

        table.add_column(Column(newcol, name=add_key))
    else:
        table[add_key]=table[add_key].astype(np.dtype('<U320'))
        for i,row in enumerate(table):
            row_data = get_data[add_key](row, table=table)
            table[i][add_key] = str(row_data)

    return(table)


def check_NED(row):
    coord = get_coord(row)
    if not coord:
        return({'error': 'no coordinate'})

    try:
        data = Ned.query_region(coord, radius=8*u.arcmin)
    except requests.exceptions.ConnectionError:
        return({'error': 'failed search'})

    # Get only data points with redshifts
    mask = (~data['Redshift'].mask) & (data['Type']==b'G')
    data = data[mask]

    if len(data)==0:
        return({'error': 'no candidates'})

    # Assume best is roughly smallest physical separation
    # Should be roughly in kpc
    best = np.argmin(data['Separation']*data['Redshift']*1262.)
    best = data[best]
    url='https://ned.ipac.caltech.edu/?q=byname&objname={0}'
    url=url.format(best['Object Name'].replace(' ','+'))
    output = {'name': best['Object Name'], 'ra': best['RA'],
              'DEC': best['DEC'], 'redshift': best['Redshift'],
              'error': None,
              'separation': best['Separation']*best['Redshift']*1262.,
              'url': url}

    return(output)

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
    result = result['GWGC','HyperLEDA','z','RAJ2000','DEJ2000','Dist','BMAG',
        'Bmag']

    for i,row in enumerate(result):
        if (result['BMAG'].mask[i] and ~result['Bmag'].mask[i] and
            ~result['Dist'].mask[i]):
            dm=5.0*np.log10(result['Dist'][i])+25.0
            result['BMAG'][i]=result['Bmag'][i]-dm
            result['BMAG'].mask[i]=False


    # Add objects because GLADE does not contain M31, M33
    result.add_row(['M31','M31',-0.001001,10.6847929,
        41.2690650,0.731,-20.10,4.22])
    result.add_row(['M33','M33',-0.000597,23.4620417,
        30.6602222,0.832,-18.31,6.29])

    mask = ~result['GWGC'].mask & ~result['z'].mask & (result['GWGC']!='---')
    result = result[mask]

    sep = []
    # Approximate separation in kpc
    distmask = result['Dist'].mask
    for i,row in enumerate(result):
        c1 = SkyCoord(row['RAJ2000'], row['DEJ2000'], unit='deg')
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
    newcol = Column(Lb*10**(
        -0.4*(result['BMAG']-5.48))/result['Separation']**2,
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

def get_hst_data(row, table=None):
    coord = get_coord(row)
    radius = 2.0 * u.arcmin

    productlist = get_productlist(coord, radius)

    if not productlist:
        return({'error': 'bad dict'})
    else:
        return({'data': productlist})

def get_osc_data(row, table=None):

    output = {}
    name = row['Name']
    if is_number(name[0]):
        tries = ['SN'+name,'AT'+name,name]
    else:
        tries = [name,'SN'+name,'AT'+name]

    for objname in tries:
        if 'ASASSN' in objname and 'ASASSN-' not in objname:
            objname = objname.replace('ASASSN','ASASSN-')

        daturl = 'https://sne.space/astrocats/astrocats/supernovae/output/json/'
        daturl += objname+'.json'

        try:
            r = requests.get(daturl)
        except requests.exceptions.ConnectionError:
            return(output)
        if r.status_code==200:
            data = json.loads(r.content)
            data = data[objname]

            if 'alias' in data.keys():
                for alias in data['alias']:
                    val=alias['value']
                    newurl='https://sne.space/astrocats/astrocats/supernovae/output/json/'
                    newurl+=val+'.json'
                    r = requests.get(newurl)
                    if r.status_code==200:
                        newdata = json.loads(r.content)
                        newdata = newdata[val]
                        if len(newdata.keys())>len(data.keys()):
                            daturl=newurl
                        data.update(newdata)


            url=daturl.replace('/astrocats/astrocats/supernovae/output/json',
                '/sne')
            url=url.replace('.json','')

            output['url']=url
            output['json']=daturl

            # Split off relevant metadata into a single dictionary
            if 'photometry' in data.keys():
                table = None
                for phot in data['photometry']:
                    photkeys=['time', 'band', 'magnitude', 'e_magnitude']
                    if all([k in phot.keys() for k in photkeys]):
                        newrow = []
                        for k in photkeys:
                            if k in phot.keys():
                                if not isinstance(phot[k], list):
                                    newrow.append(phot[k])
                        if len(newrow)==len(photkeys):
                            if not table:
                                newrow=[[el] for el in newrow]
                                table = Table(newrow, names=photkeys)
                            else:
                                table.add_row(newrow)

                output['photometry'] = table

            if 'discoverdate' in data.keys():
                t = Time(data['discoverdate'][0]['value'].replace('/','-'))
                date_str = t.datetime.strftime('%Y-%m-%d %H:%M:%S')
                output['discovery_date'] = date_str

            if 'redshift' in data.keys():
                output['redshift'] = data['redshift'][0]['value']

            if 'spectra' in data.keys():
                output['spectra'] = data['spectra']

            if 'sources' in data.keys():
                sources = {}
                for source in data['sources']:
                    if 'bibcode' in source.keys() and 'alias' in source.keys():
                        sources[source['alias']]=source['bibcode']

                output['sources'] = sources

            if 'claimedtype' in data.keys():
                typs=[]
                classes=[]
                for typ in data['claimedtype']:
                    typs.append(typ['source'])
                    classes.append(typ['value'])

                output['classification_reference']=','.join(typs)
                output['classification']=','.join(classes)

            if 'discoverdate' in data.keys():
                disc_refs=[]
                for dat in data['discoverdate']:
                    disc_refs.append(dat['source'])

                output['discovery_reference']=','.join(disc_refs)

            break

    return(output)

def get_classification(row, table=None):

    tns_class = check_dict(table.meta, [row['Name'],'tns','object_type','name'])
    osc_class = check_dict(table.meta, [row['Name'],'osc','classification'])

    objtype=''
    if osc_class:
        osc_class = osc_class.split(',')[0]
        osc_class = osc_class.strip()
        objtype = osc_class

    if tns_class and (not objtype or objtype.upper()=='CANDIDATE'):
        tns_class = tns_class.strip()
        objtype = tns_class

    if objtype:
        if objtype.upper() in ['NOVA','LBV','LRN','CANDIDATE','ILRT',
            'LBV TO IIN','TDE']:
            return(objtype)

        if not objtype.upper().startswith('SN'):
            objtype = 'SN '+objtype
        return(objtype)

    return('')

def get_classification_mask(table):

    if 'Classification' not in table.keys():
        return(None)

    class_mask = []
    for row in table:
        if row['Classification'].upper() in ['CANDIDATE','NOVA','SN','']:
            class_mask.append(True)
        elif not row['Classification'].strip():
            class_mask.append(True)
        else:
            class_mask.append(False)
    class_mask = np.array(class_mask)

    pre_exp_mask = []
    for row in table:
        disc_time = Time(row['Discovery Date'])
        data = check_dict(table.meta, [row['Name'],'hst','data'])
        if data:
            if 'start_time' not in row.colnames:
                pre_exp_mask.append(True)
                continue
            times = Time(row['start_time'], format='mjd')
            mask = times < disc_time
            if len(data[mask])==0:
                pre_exp_mask.append(True)
            else:
                pre_exp_mask.append(False)
        else:
            pre_exp_mask.append(True)
    pre_exp_mask = np.array(pre_exp_mask)

    mask = class_mask & pre_exp_mask

    return(mask)


def get_ads_data(row, table=None):

    baseurl='https://api.adsabs.harvard.edu/v1/search/query?q={0}'
    baseurl+='&fq=database:astronomy&fl=bibcode,title'

    if is_number(row['Name'][0]):
        tries = [row['Name'],'AT'+row['Name'],'SN'+row['Name']]
    else:
        tries = [row['Name']]

    return_val = {}
    for t in tries:
        url = baseurl.format(row['Name'])
        try:
            r = requests.get(url, headers=ads_headers)
        except requests.exceptions.ConnectionError:
            return({'error': 'timeout error'})

        if r.status_code==200 and 'response' in r.json().keys():
            output = r.json()['response']
            if 'docs' in output.keys() and  output['numFound']>0:
                if not return_val:
                    return_val=output
                else:
                    bibcodes = [v['bibcode'] for v in return_val['docs']]
                    for val in output['docs']:
                        if val['bibcode'] not in bibcodes:
                            return_val['docs'].append(val)
                            return_val['numFound']+=1

            return_val['error']=None

    return(return_val)

def get_host_distance(row, table=None):

    name = check_dict(table.meta, [row['Name'],'ned','name'])

    df = []
    name = row['Host']

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

    return({})

def get_best_distance(row, table=None):

    table = check_dict(table.meta, [row['Name'],'distance','table'])

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


def get_host(row, table=None):

    name = check_dict(table.meta, [row['Name'], 'ned', 'name'])
    if name: return(name)

    return('')

def format_date(date):
    t = Time(date)
    date_str = t.datetime.strftime('%Y-%m-%d %H:%M')
    sec_str  = str(int(np.round(float(t.datetime.strftime('%S'))))).zfill(2)
    date_str = date_str+':'+sec_str

    return(date_str)


def get_disc_date(row, table=None):

    date = check_dict(table.meta, [row['Name'], 'tns','discoverydate'])
    if (date and date.strip() and ('20' in date or '19' in date)):
        return(format_date(date))

    date = check_dict(table.meta, [row['Name'], 'osc', 'discovery_date'])
    if (date and date.strip() and ('20' in date or '19' in date)):
        return(format_date(date))

    date = check_dict(table.meta, [row['Name'], 'yse', 'discoverydate'])
    if (date and date.strip() and ('20' in date or '19' in date)):
        return(format_date(date))

    return('')

def get_redshift(row, table=None):

    ned_redshift = check_dict(table.meta, [row['Name'], 'ned', 'redshift'])
    if ned_redshift and is_number(ned_redshift):
        return(round_to_n(ned_redshift, 6))

    osc_redshift = check_dict(table.meta, [row['Name'], 'osc', 'redshift'])
    if osc_redshift and is_number(osc_redshift):
        return(round_to_n(osc_redshift, 6))

    return('')

def get_yse_targets():
    url = 'https://ziggy.ucolick.org/yse/explorer/160/download?format=csv'
    try:
        r = requests.get(url)
    except requests.exceptions.ConnectionError:
        return(None)
    if r.status_code==200:
        table = ascii.read(r.text)
        for key in table.keys():
            if 'name' in key:
                table.rename_column(key, 'name')
        return(table)
    else:
        return(None)

def get_yse_target_photometry():
    url = 'https://ziggy.ucolick.org/yse/explorer/218/download?format=csv'
    try:
        r = requests.get(url)
    except requests.exceptions.ConnectionError:
        return(None)
    if r.status_code==200:
        table = ascii.read(r.text)
        for key in table.keys():
            if 'name' in key:
                table.rename_column(key, 'name')
        return(table)
    else:
        return(None)

def add_yse_targets(sndata):

    yse_targets = get_yse_targets()
    if not yse_targets: return(sndata)

    for row in yse_targets:

        if any([row['name'] in sndata[key]['Name'] for key in sndata.keys()]):
            continue

        specclass = row['TNS_spec_class']
        add_key = 'Other'
        if specclass=='SN Ia': add_key = 'Type Ia'
        if specclass=='SN IIb': add_key = 'Type IIb'
        if specclass=='SN Ib': add_key = 'Type Ib/c'
        if specclass=='SN Ic': add_key = 'Type Ib/c'
        if specclass=='SN Ib Pec': add_key = 'Type Ib/c'
        if specclass=='SN Ic-BL': add_key = 'Type Ib/c'
        if specclass=='SN IIP': add_key = 'Type II-P/II-L'
        if specclass=='SN IIL': add_key = 'Type II-P/II-L'
        if specclass=='SN II': add_key = 'Type II-P/II-L'
        if specclass=='SN IIn': add_key = 'Type IIn'
        if specclass=='LBV': add_key = 'Type IIn'

        coord = SkyCoord(row['transient_RA'], row['transient_Dec'], unit='deg')
        ra, dec = coord.to_string(style='hmsdms', sep=':', precision=3).split()
        date_str = Time(row['disc_date']).datetime.strftime('%Y-%m-%d %H:%M:%S')

        add_table = sndata[add_key]
        new_row = []
        for tab_key in add_table.keys():
            if tab_key=='Name': new_row.append(row['name'])
            elif tab_key=='RA': new_row.append(ra)
            elif tab_key=='Dec': new_row.append(dec)
            elif tab_key=='Redshift': new_row.append(str(row['redshift']))
            elif tab_key=='Discovery Date': new_row.append(date_str)
            else: new_row.append('')

        add_table.add_row(new_row)
        sndata[add_key]=add_table

        # If yes, send an alert
        from_addr = email_args['from_addr']
        smtpserver = email_args['smtpserver']
        gmail_login = email_args['gmail_login']
        gmail_password = email_args['gmail_password']
        email = email_args['to_addr']

        email_summary = 'New Supernova Progenitor Target<br>'
        email_summary += '{0}<br>'.format(row['name'])
        email_summary += '{0}<br>'.format(add_key)
        email_summary += 'z={0}<br>'.format(row['redshift'])
        email_summary += '{0}<br>'.format(date_str)

        resp=sendEmail(from_addr, email, email_args['subject'],
            email_summary, gmail_login, gmail_password, smtpserver)

    return(sndata)

def add_metadata(table, method, redo=False):

    for row in table:
        method=method.lower()

        print('Getting {0} data for {1}'.format(method, row['Name']))

        if method not in get_metadata.keys():
            warning = 'WARNING: bad metadata method {0}'
            print(warning.format(method))
            continue

        if row['Name'] not in table.meta.keys():
            redo = True
            print('No name')
        elif method not in table.meta[row['Name']].keys():
            redo = True
            print(f'{method}, No method')
        elif (table.meta[row['Name']][method] is None or
              table.meta[row['Name']][method]=={}):
            redo = True
            print('Method is None')
        elif ('error' in table.meta[row['Name']][method].keys() and
            (table.meta[row['Name']][method]['error'] in ['timeout',
                'timeout data', 403, 'bad dict','no candidates'])):
            redo = True
            print('Error',table.meta[row['Name']][method]['error'])
        else:
            print(table.meta[row['Name']][method].keys())

        if redo:
            print('Redoing (or doing, whatever)...')
            metadata = get_metadata[method](row, table=table)
            if row['Name'] not in table.meta.keys():
                table.meta[row['Name']]={method: metadata}
            else:
                table.meta[row['Name']][method]=metadata

    return(table)

get_metadata = {'osc': get_osc_data,
                'tns': get_tns_data,
                'yse': get_yse_data,
                'ned': check_glade,
                'distance': get_host_distance,
                'ads': get_ads_data,
                'hst': get_hst_data,}

def get_ads_ref(row, table, typ):

    ads_dict = check_dict(table.meta, [row['Name'], 'ads'])
    if not ads_dict:
        return('')

    if 'numFound' not in ads_dict.keys() or ads_dict['numFound']==0:
        return('')

    url = 'https://ui.adsabs.harvard.edu/abs/{0}'

    refs = ads_dict['docs']
    idx = []
    if typ=='discovery':
        idx = [i for i,val in enumerate(refs)
            if ((len(val['title'])>0) and
                ('discovery' in val['title'][0].lower()))]
    elif typ=='classification':
        idx = [i for i,val in enumerate(refs)
            if ((len(val['title'])>0) and
                ('classification' in val['title'][0].lower()) or
                ('spectroscopic' in val['title'][0].lower()) or
                ('spectroscopy' in val['title'][0].lower()))]

    if len(idx)>0:
        idx = np.max(idx)
        url = url.format(refs[idx]['bibcode'])
        return(url)

    return('')

def get_distance(row, table=None):

    out_fmt = '{0} ({1})'
    output = get_best_distance(row, table=table)
    if 'distance' in output.keys() and 'distance_error' in output.keys():
        out = out_fmt.format(output['distance'], output['distance_error'])
        return(out)

    return('')

def get_host_url(row, table=None):
    url = check_dict(table.meta, [row['Name'], 'ned', 'url'])
    if url: return(url)

    return('')

def get_distance_method(row, table=None):
    output = get_best_distance(row, table=table)
    if output and 'method' in output.keys():
        return(output['method'])

    return('')

def get_distance_ref(row, table=None):
    output = get_best_distance(row, table=table)
    if output and 'reference' in output.keys():
        ref = output['reference']
        return('https://ui.adsabs.harvard.edu/abs/{0}'.format(ref))

    return('')

def get_osc_source(row, table, refname):
    sources = check_dict(table.meta, [row['Name'],'osc','sources'])
    ref = check_dict(table.meta, [row['Name'],'osc',refname])

    if sources and ref:
        # Resolve references into bibcodes
        all_sources = []
        for r in ref.split(','):
            if r in sources.keys():
                all_sources.append(sources[r])

        if len(all_sources)>0:
            all_sources = sorted(all_sources, key=lambda s: s[0:4])
            # Take the most recent source
            source = all_sources[-1]

            url='https://ui.adsabs.harvard.edu/abs/{0}'
            return(url.format(source))

    return('')

def get_tns_source(row, tns_cert_type):

    if tns_cert_type.lower()=='discovery':
        url='https://www.wis-tns.org/object/{0}/discovery-cert'
    elif tns_cert_type.lower()=='classification':
        url='https://www.wis-tns.org/object/{0}/classification-cert'
    else:
        return('')

    source = ''

    url = url.format(row['Name'])
    headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) '+\
        'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 '+\
        'Safari/537.36'}
    try:
        r = requests.get(url, headers=headers)
    except requests.exceptions.ConnectionError:
        return(source)

    if r.status_code==200:
        reg = '\"https://ui\.adsabs\.harvard\.edu/search/q=.*\"'
        matches = re.findall(reg, r.text)
        if len(matches)>0:
            source = matches[0].split()[0].replace('\"','')

    return(source)

def get_discovery_ref(row, table=None):
    osc_source = get_osc_source(row, table, 'discovery_reference')
    if table:
        tns_source = check_dict(table.meta, [row['Name'],'tns',
            'discovery_reference'])
    ads_source = get_ads_ref(row, table, 'discovery')

    if osc_source: return(osc_source)
    if tns_source: return(tns_source)
    if ads_source: return(ads_source)
    return('')

def get_type_ref(row, table=None):
    osc_source = get_osc_source(row, table, 'classification_reference')
    if table:
        tns_source = check_dict(table.meta, [row['Name'],'tns',
            'classification_reference'])
    ads_source = get_ads_ref(row, table, 'classification')

    if osc_source: return(osc_source)
    if tns_source: return(tns_source)
    if ads_source: return(ads_source)
    return('')

def get_spectra(row, table):

    all_spectra = []

    spec_osc = check_dict(table.meta, [row['Name'],'osc','spectra'])
    if spec_osc:
        for spec in spec_osc:
            if ('u_wavelengths' in spec.keys() and
                'time' in spec.keys() and
                'u_time' in spec.keys() and
                spec['u_time']=='MJD' and
                spec['u_wavelengths'].lower()=='angstrom'):
                time = Time(float(spec['time']), format='mjd')
                data = spec['data']

                spectrum = [[float(d[0]),float(d[1])] for d in data]
                all_spectra.append({'time': time, 'data': spectrum})

    return(all_spectra)

def gather_type(sndata, sn_type):

    type_index = {'Type Ia': ['SN Ia', 'SN Ia Pec', 'SN Ia-91bg',
                    'SN Ia-02cx','SN Ia-91bg-like','SN Ia-91T',
                    'SN Ia-pec','SN I','SN Iax[02cx-like]'],
                  'Type Ib/c': ['SN Ib','SN Ic','SN Ib/c','SN Ic BL',
                    'SN Ic-BL','SN Ib-pec','SN Ibn'],
                  'Type II-P/II-L': ['SN II','SN II P','SN IIP'],
                  'Type IIn': ['LBV','SN IIn','SN IIn/LBV','SN IIn-pec/LBV',
                    'SN LBV to IIn','LBV to IIn'],
                  'Type IIb': ['SN IIb'],
                  'Other': ['ILRT','LRN','SN','Nova','TDE'],}

    new_table = sndata[sn_type].copy()[:0]
    new_table.meta = {}

    acceptable_types = type_index[sn_type]+['Candidate','']
    for row in sndata[sn_type]:
        if row['Classification'] in acceptable_types:
            if row['Name'] not in new_table['Name']:
                new_table.add_row(row)
                new_table.meta[row['Name']]=sndata[sn_type].meta[row['Name']]

    acceptable_types = type_index[sn_type]
    for typ in sndata.keys():
        if typ==sn_type: continue
        for row in sndata[typ]:
            if row['Classification'] in acceptable_types:
                if row['Name'] not in new_table['Name']:
                    new_table.add_row(row)
                    new_table.meta[row['Name']]=sndata[typ].meta[row['Name']]

    return(new_table)

def get_hst(row, table=None):
    data = check_dict(table.meta, [row['Name'],'hst','data'])
    if data:
        exptime = np.sum(data['exptime'].data)
        exptime = '%10.2f'%exptime
        exptime = exptime.strip()
        return(exptime)

    return(0)

get_data = {'OSC': get_osc,
            'TNS': get_tns,
            'YSEPZ': get_yse,
            'Host': get_host,
            'Discovery Date': get_disc_date,
            'Distance': get_distance,
            'Redshift': get_redshift,
            'HST': get_hst,
            'NED': get_host_url,
            'Distance Method': get_distance_method,
            'Ref. (Distance)': get_distance_ref,
            'Ref. (Discovery)': get_discovery_ref,
            'Ref. (Classification)': get_type_ref,
            'Classification': get_classification,}

def get_productlist(coord, search_radius):

    productlist = None

    # Check for coordinate and exit if it does not exist
    if not coord:
        error = 'ERROR: coordinate was not provided.'
        return(productlist)

    try:
        obsTable = Observations.query_region(coord, radius=search_radius)
    except (astroquery.exceptions.RemoteServiceError,
        requests.exceptions.ConnectionError,
        astroquery.exceptions.TimeoutError,
        requests.exceptions.ChunkedEncodingError):
        error = 'ERROR: MAST is not working currently working\n'
        error += 'Try again later...'
        print(error)
        return(productlist)

    # Get rid of all masked rows (they aren't HST data anyway)
    obsTable = obsTable.filled()

    # Construct masks for telescope, data type, detector, and data rights
    masks = []
    masks.append([t.upper()=='HST' for t in obsTable['obs_collection']])
    masks.append([p.upper()=='IMAGE' for p in obsTable['dataproduct_type']])
    masks.append([any(l) for l in list(map(list,zip(*[[det in inst.upper()
                for inst in obsTable['instrument_name']]
                for det in ['ACS','WFC','WFPC2']])))])

    # Get rid of short exposures (defined as 15s or less)
    masks.append([t > 15. for t in obsTable['t_exptime']])

    # Apply the masks to the observation table
    mask = [all(l) for l in list(map(list, zip(*masks)))]
    obsTable = obsTable[mask]

    # Iterate through each observation and download the correct product
    # depending on the filename and instrument/detector of the observation
    for obs in obsTable:
        try:
            productList = Observations.get_product_list(obs)
            # Ignore the 'C' type products
            mask = productList['type']=='S'
            productList = productList[mask]
        except:
            error = 'ERROR: MAST is not working currently working\n'
            error += 'Try again later...'
            print(error)
            return(productlist)

        instrument = obs['instrument_name']
        s_ra = obs['s_ra']
        s_dec = obs['s_dec']
        exptime = obs['t_exptime']
        filt = obs['filters']
        start = obs['t_min']

        instcol = Column([instrument]*len(productList), name='instrument_name')
        racol = Column([s_ra]*len(productList), name='ra')
        deccol = Column([s_dec]*len(productList), name='dec')
        expcol = Column([exptime]*len(productList), name='exptime')
        filtcol = Column([filt]*len(productList), name='filter')
        startcol = Column([start]*len(productList), name='start_time')

        productList.add_column(instcol)
        productList.add_column(racol)
        productList.add_column(deccol)
        productList.add_column(expcol)
        productList.add_column(filtcol)
        productList.add_column(startcol)

        for prod in productList:
            filename = prod['productFilename']

            if (('c0m.fits' in filename and 'WFPC2' in instrument) or
                ('c1m.fits' in filename and 'WFPC2' in instrument) or
                ('c0m.fits' in filename and 'PC/WFC' in instrument) or
                ('c1m.fits' in filename and 'PC/WFC' in instrument) or
                ('flc.fits' in filename and 'ACS/WFC' in instrument) or
                ('flt.fits' in filename and 'ACS/HRC' in instrument) or
                ('flc.fits' in filename and 'WFC3/UVIS' in instrument) or
                ('flt.fits' in filename and 'WFC3/IR' in instrument)):

                if not productlist:
                    productlist = Table(prod)
                else:
                    productlist.add_row(prod)

    if not productlist:
        return(productlist)

    downloadFilenames = []
    for prod in productlist:
        filename = prod['productFilename']

        # Cut down new HST filenames that start with hst_PROGID
        filename = '_'.join(filename.split('_')[-2:])
        downloadFilenames.append(filename)

    productlist.add_column(Column(downloadFilenames, name='downloadFilename'))

    # Check that all files to download are unique
    if productlist and len(productlist)>1:
        productlist = unique(productlist, keys='downloadFilename')

    # Sort by obsID in case we need to reference
    productlist.sort('obsID')

    return(productlist)
