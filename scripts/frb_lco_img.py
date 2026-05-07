import util
import sheetproc
import sys
import datetime
import numpy as np
import math
from scipy import signal, integrate
from lcogt import lcogt
from astropy.time import Time, TimeDelta
from astropy.table import unique, vstack, Table
from astropy.coordinates import SkyCoord
from astropy import units as u

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from scipy.interpolate import interp1d

chime = EarthLocation(lat=49.3211111111*u.deg, lon=-119.623888889*u.deg, height=0*u.m)

from astropy.coordinates import get_sun
def get_sun_altitudes(obs, t=Time('2022-02-01 00:00:00'),
    tend=Time('2022-08-01 00:00:00'), tdelta=TimeDelta(300*u.s)):

    N=int((tend-t)/tdelta)
    deltas = np.linspace(0, N, N)*tdelta
    alltimes=t+deltas-TimeDelta(3600*u.s)

    M=int(86400*u.s/tdelta)
    n=int(N/M)

    sun_times = get_sun(alltimes)
    sun_altitudes = sun_times.transform_to(AltAz(location=obs)).alt.degree
    sun_altitudes = np.reshape(sun_altitudes, (n, M))

    return(sun_altitudes)

def get_altitudes(coord, obs, t=Time('2022-02-01 00:00:00'),
    tend=Time('2022-08-01 00:00:00'), tdelta=TimeDelta(300*u.s)):

    N=int((tend-t)/tdelta)
    deltas = np.linspace(0, N, N)*tdelta
    alltimes=t+deltas

    M=int(86400*u.s/tdelta)
    n=int(N/M)

    altitudes=coord.transform_to(AltAz(obstime=alltimes,location=obs)).alt.degree
    altitudes=np.reshape(altitudes, (n, M))

    return(altitudes)

# Analysis functions
# Given a set of MJDs, calculate the period and t0 value for sinusoid
# as well as periodogram frequencies and power spectrum
def get_period_t0(mjds):
    ran = np.linspace(np.min(mjds), np.max(mjds), 1000)
    sig = np.zeros(len(ran))
    scale = ran[1]-ran[0]

    # Rebin data into uniform distribution of dates
    for i,el in enumerate(ran):
        if i < len(ran)-1:
            mask = (mjds > el) & (mjds < ran[i+1])
            sig[i] = len(mjds[mask])

    f, Pxx_den = signal.periodogram(sig)
    # Rescale frequency to match date distribution
    f = f/scale

    # Get rid of negative frequency values
    mask = f > 0
    f = f[mask] ; Pxx_den = Pxx_den[mask]

    # Normalize power spectrum
    Pxx_den = Pxx_den / np.max(Pxx_den)

    # Get period for peak of periodogram
    period = 1./f[np.argmax(Pxx_den)]

    # To get t0 of sinusoid, determine t0 in distribution from 0->period
    # where sum of dates at that time in sinusoid is maximized
    t_init = np.linspace(0, period, 1000)
    fitting = np.zeros(len(t_init))
    for i,t in enumerate(t_init):
        fitting[i] = np.sum([math.cos(2*np.pi*(Time(d).mjd-t)/period)
            for d in dates])

    t0 = t_init[np.argmax(fitting)]

    return(period,t0,f,Pxx_den)

band_priority=['r','rp','r-ZTF','R','orange','V','i','ip','I','g','gp','g-ZTF',
    'cyan','B']
lco = lcogt('/home/ckilpatrick/scripts/shibboleth')
sdate = Time(datetime.datetime.now()) - TimeDelta(28, format='jd')
rg_img = lco.get_requestgroups(propid=lco.proposals,
    itype='1M0-SCICAM-SINISTRO', sdate=sdate)

frb_targets = Table([['FRB180916','FRB200120E','FRB181030A'],
    ['01:58:00','09:57:56.7','10:34:51.2'],
    ['+65:43:00','+68:49:32','+73:44:38']],names=('name','ra','dec'))

table = Table.read('frb180916.dat', format='ascii',
    names=('name','date','time','dm','dm_err','snr'))
dates = [datetime.datetime.strptime(d+'T'+t, "%Y-%m-%dT%H:%M:%S.%f")
    for d,t in zip(table['date'],table['time'])]
mjds = np.array([Time(d).mjd for d in dates])
ago  = Time(datetime.datetime.utcnow()) - TimeDelta(720, format='jd')
mask = mjds > ago.mjd
mjds = mjds[mask]
period, t_rel, f, Pxx_den = get_period_t0(mjds)

print(period, t_rel)

for i in np.arange(10):
    now = Time(datetime.datetime.utcnow()).mjd
    periods = np.floor((now-t_rel) / period)

    mjd = (periods-1+i)*period + t_rel
    t = Time(mjd, format='mjd')
    print(t.datetime.strftime('%Y-%m-%d %H:%M:%S'))

def main(redo=False):

    for target in frb_targets:

        targcoord = SkyCoord(target['ra'], target['dec'], unit=(u.hour,u.deg))
        name = target['name']

        ra = targcoord.ra.degree
        dec = targcoord.dec.degree

        now = Time(datetime.datetime.utcnow())
        interval = TimeDelta(1, format='jd')
        tend = now + interval

        tdelta = TimeDelta(5*u.s)
        chime_altitudes = get_altitudes(targcoord, chime, t=now, tend=tend,
                tdelta=tdelta)
        delta_t = np.argmax(chime_altitudes[0])*tdelta
        best_time = now + delta_t

        print(f'Best time to observe {name} with CHIME: {best_time}')

        start_time = best_time - TimeDelta(8*60, format='sec')

        needs_obs = True
        cadence = 9.5
        if target['name']=='FRB180916':
            now = Time(datetime.datetime.utcnow()).mjd
            periods = np.floor((now-t_rel) / period)
            rel = np.abs(now - t_rel - periods * period)
            if rel < 2.5:
                cadence = 0.75
            else:
                cadence = 99.0

        print(f'Cadence is: {cadence}')

        # Get LCO data to determine if we need a new observation
        for request in rg_img:
            # If PENDING observation, don't need a new one
            if request['state']=='PENDING':
                window = request['requests'][0]['windows'][0]
                start = Time(window['start']).mjd
                end = Time(window['end']).mjd
                now = Time(datetime.datetime.utcnow()).mjd
                if start < now and end < now:
                    continue
                elif start > now+cadence and end > now+cadence:
                    continue
                for config in request['requests'][0]['configurations']:
                    if ('ra' not in config['target'].keys() or
                        'dec' not in config['target'].keys()):
                        continue
                    coord = SkyCoord(config['target']['ra'],
                        config['target']['dec'], unit='deg')
                    if coord.separation(targcoord).degree < 0.3:
                        needs_obs = False

            # Use start time as a proxy for observation time and cut on observations
            # that are older than the now - cadence
            time = Time(request['requests'][0]['modified']).mjd
            now = Time(datetime.datetime.now()).mjd
            if now - time > cadence:
                continue

            # For records more recent than cadence, check COMPLETED and PENDING
            # observations.
            if request['state']=='COMPLETED':
                window = request['requests'][0]['windows'][0]
                start = Time(window['start']).mjd
                end = Time(window['end']).mjd
                now = Time(datetime.datetime.utcnow()).mjd
                if start < now-cadence and end < now-cadence:
                    continue
                for config in request['requests'][0]['configurations']:
                    coord = SkyCoord(config['target']['ra'],
                        config['target']['dec'], unit='deg')
                    # If the separation between the recent or pending observation
                    # and the target coordinates is smaller than threshold, then
                    # record that this target does not need to be re-observed
                    if coord.separation(targcoord).degree < 0.3:
                        needs_obs = False

        if needs_obs:
            ra = targcoord.ra.degree
            dec = targcoord.dec.degree

            now = Time(datetime.datetime.utcnow())
            interval = TimeDelta(1, format='jd')
            tend = now + interval

            tdelta = TimeDelta(5*u.s)
            chime_altitudes = get_altitudes(targcoord, chime, t=now, tend=tend,
                tdelta=tdelta)
            delta_t = np.argmax(chime_altitudes[0])*tdelta
            best_time = now + delta_t

            start_time = best_time - TimeDelta(8*60, format='sec')

            response = lco.make_obs_request(target['name'], ra, dec, 18.0,
                propidx=0, strategy = 'photometry-frb-time-critical',
                start=start_time.datetime)

            if response and 'requests' in response.keys():
                if 'non_field_errors' in response['requests'][0].keys():
                    message = '{target} cannot be scheduled due to availability.'
                    print(response['requests'][0])
                    print(message.format(target=target['name']))
                else:
                    message = 'Successfully scheduled {target}'
                    print(message.format(target=target['name']))
            else:
                message = '{target} exposures not possible with current settings'
                print(message.format(target=target['name']))
                print(response)
        else:
            message = '{target} does not need a new observation.'
            print(message.format(target=target['name']))

main(redo=False)

