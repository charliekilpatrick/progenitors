from astropy.io import ascii
import numpy as np
from scipy import interpolate
from scipy.integrate import simps
import pysynphot as S
import sys
import pickle

dustdir = 'data/dust/'

# RSG models for range of temperatures
temp=[2600,2800,3000,3200,3300,3400,3500,3600,3700,3800,3900,4000,
            4250,4500,5000,6000,7000,8000]
temp = np.array(temp)

# This is the total integrated flux of a source with Lum = 1 Lsun in units of
# erg/s/cm2.  Used to renormalize the RSG spectrum.  This is equivalent to
# 3.839e33/(4 * pi * (10 * 3.08568025e18)^2) for 1 Lsol at 10 pc.
FLUX_SCALE = 3.22398177e-7

# Dust-to-gas mass ratio assumption
DUST_TO_GAS = 0.01

data10 = np.loadtxt(dustdir+'data10.dat', unpack=True, dtype=float)
#data05 = np.loadtxt(dustdir+'data05.dat', unpack=True, dtype=float)
#data0_5 = np.loadtxt(dustdir+'data0-5.dat', unpack=True, dtype=float)
#data00 = np.loadtxt(dustdir+'data00.dat', unpack=True, dtype=float)

kappa = np.loadtxt(dustdir+'dust01_trans.dat')
wavelength = np.loadtxt(dustdir+'wavelength.dat')

def rebin(a, newshape):
    newarray = np.zeros(newshape)
    curr = 0
    for i in np.arange(newshape):
        s = slice(curr, curr+int(a.shape[0]/newshape), 1)
        newarray[i] = np.mean(a[s])
        curr += int(a.shape[0]/newshape)
    return(newarray)

wavelength = rebin(wavelength, 7748)
rsg_10 = interpolate.interp2d(wavelength, temp, data10)
#rsg_05 = interpolate.interp2d(wavelength, temp[:data05.shape[0]], data05)
#rsg_0_5 = interpolate.interp2d(wavelength, temp[:data0_5.shape[0]], data0_5)
#rsg_00 = interpolate.interp2d(wavelength, temp[:data00.shape[0]], data00)

# kappa_V, the opacity in V-band calculated from dust01_trans.dat
bp = S.ObsBandpass('johnson,V')
sp1 = S.ArraySpectrum(wavelength, kappa, 
    fluxunits='counts', waveunits='angstrom')
sp2 = S.ArraySpectrum(wavelength, np.array([1.0]*len(wavelength)), 
    fluxunits='counts', waveunits='angstrom')
obs1 = S.Observation(sp1, bp, binset=wavelength)
obs2 = S.Observation(sp2, bp, binset=wavelength)
kappa_V=obs1.effstim('counts')/obs2.effstim('counts')

# Normalize kappa for calculation below so we can derive correct luminosity
kappa = kappa/np.max(kappa)

# Graphite, Rout/Rin=2
def get_avg2(x, p):
    l=x*1.0e-4
    t=p/10.0
    tgra1 = (0.500446 + 1.795729*t - 1.877658*t*t + 0.852820*t**3 - 0.141635*t**4)
    tgra2 = (4.318269 - 14.236698*t + 13.804110*t*t - 5.991369*t**3+0.959539*t**4)/l
    tgra3 = (-5.114167 + 32.462564*t - 26.895305*t*t + 10.197398*t**3 - 1.414338*t**4)/l**2
    tgra4 = (3.384105 - 21.107633*t + 12.229167*t*t - 2.172318*t**3 - 0.149866*t**4)/l**3
    tgra5 = (-1.059677 + 5.553703*t - 1.527415*t*t - 0.881450*t**3+0.391763*t**4)/l**4
    tgra6 = (0.121772 - 0.518968*t - 0.048566*t*t + 0.248002*t**3 - 0.077054*t**4)/l**5
    agraphite2 = (tgra1 + tgra2 + tgra3 + tgra4 + tgra5 + tgra6)*t*l**(-1.375229)

    return(agraphite2)

def get_gvg2(x, p):
    l=x*1.0e-4
    s=(p[0]/10.0)**(1/2)
    ggra1 = (0.091091-2.031641*s + 5.049221*s*s-6.399577*s**3+2.086527*s**4)
    ggra2 = (-0.705931 + 3.566326*s-40.821335*s*s + 48.029644*s**3-14.815801*s**4)/l
    ggra3 = (1.604726-22.761479*s + 90.223336*s*s-88.179577*s**3 + 25.041643*s**4)/l**2
    ggra4 = (-1.240894 + 20.383229*s-63.820599*s*s + 57.434801*s**3-15.475733*s**4)/l**3
    ggra5 = (0.392600-6.856117*s + 19.106438*s*s-16.340312*s**3+4.240107*s**4)/l**4
    ggra6 = (-0.044443 + 0.800158*s-2.078890*s*s + 1.716855*s**3-0.433098*s**4)/l**5
    ggraphite2 = -(ggra1 + ggra2 + ggra3 + ggra4 + ggra5 + ggra6)*s/(1 + l**(3.608885))

    return(ggraphite2)

# Graphite, Rout/Rin=10
def get_avg10(x, p):
    l=x*1.0e-4
    t=p[0]/10.0
    tgra1 = (0.760499 + 0.879164*t - 0.350748*t*t - 0.039612*t**3+0.034161*t**4)
    tgra2 = (4.061343 - 7.166933*t + 2.791544*t*t + 0.214647*t**3 - 0.233685*t**4)/l
    tgra3 = (-5.133851 + 16.344656*t - 4.283100*t*t - 1.764900*t**3+0.780217*t**4)/l**2
    tgra4 = (3.387184 - 10.066016*t - 1.260999*t*t + 4.103272*t**3 - 1.160204*t**4)/l**3
    tgra5 = (-1.052057 + 2.479576*t + 1.618868*t*t - 2.030708*t**3+0.516503*t**4)/l**4
    tgra6 = (0.120327 - 0.214118*t - 0.293914*t*t + 0.295687*t**3 - 0.071840*t**4)/l**5
    agraphite10 = (tgra1 + tgra2 + tgra3 + tgra4 + tgra5 + tgra6)*t*l**(-1.475236)

    return(agraphite10)

# Silicate, Rout/Rin=2
def get_avs2(x, p):
    l=x*1.0e-4
    t=p[0]/10.0
    tsil1 = (0.437549 - 0.446323*t + 0.648423*t*t - 0.321970*t**3+0.055555*t**4)
    tsil2 = (-0.486741 + 4.034854*t - 5.530127*t*t + 2.711095*t**3 - 0.469112*t**4)/l
    tsil3 = (1.166512 - 12.015845*t + 14.917191*t*t - 7.058630*t**3+1.204954*t**4)/l**2
    tsil4 = (-0.655682 + 13.849514*t - 14.342367*t*t + 6.165157*t**3 - 0.984406*t**4)/l**3
    tsil5 = (0.169689 - 4.956815*t + 4.137525*t*t - 1.419899*t**3+0.176166*t**4)/l**4
    tsil6 = (-0.016829 + 0.582619*t - 0.381166*t*t + 0.083593*t**3 - 0.002153*t**4)/l**5
    asilicate2 = (tsil1 + tsil2 + tsil3 + tsil4 + tsil5 + tsil6)*t*l**(-0.642318)

    return(asilicate2)

# Silicate, Rout/Rin=10
def get_avs10(x, p):
    l=x*1.0e-4
    t=p[0]/10.0
    tsil1 = (0.197398 - 0.293417*t + 0.192686*t*t - 0.041375*t**3+0.000902*t**4)
    tsil2 = (0.093593 + 2.491030*t - 1.453387*t*t + 0.239280*t**3+0.013273*t**4)/l
    tsil3 = (0.357331 - 6.883382*t + 3.239407*t*t - 0.198182*t**3 - 0.121949*t**4)/l**2
    tsil4 = (0.022567 + 7.169214*t - 1.884278*t*t - 0.728088*t**3+0.309664*t**4)/l**3
    tsil5 = (-0.065599 - 2.173130*t - 0.305525*t*t + 0.839291*t**3 - 0.223389*t**4)/l**4
    tsil6 = (0.012188 + 0.216324*t + 0.133118*t*t - 0.153598*t**3+0.036469*t**4)/l**5
    asilicate10 = (tsil1 + tsil2 + tsil3 + tsil4 + tsil5 + tsil6)*t*l**(-0.323043)

    return(asilicate10)

def get_dust(x, p, model='g2'):

    if model=='g2': dust = get_avg2(x, p)
    elif model=='s2': dust = get_avs2(x, p)
    elif model=='g10': dust = get_avg10(x, p)
    elif model=='s10': dust = get_avs10(x, p)

    dust[np.where(dust < 0)] = 0

    return(dust)

def get_rsg(scale, temp, model='10'):

    if model=='10': flux=rsg_10(wavelength, temp)
    elif model=='05': flux=rsg_05(wavelength, temp)
    elif model=='0_5': flux=rsg_0_5(wavelength, temp)
    elif model=='00': flux=rsg_00(wavelength, temp)

    normalize = simps(flux, wavelength)
    # Renormalize the RSG spectrum so it's in units of erg/s/cm2/angstrom for
    # pysynphot to interpret
    flux = scale * FLUX_SCALE * flux/normalize

    return(flux)

def get_bb_lum(p, rsg_model='10', dust_model='g2', sptype='all', masked=True):

    w = wavelength
    flux = get_rsg(p[1], p[2], model=rsg_model)

    if masked and len(w)!=5157:
        mask = (wavelength > 3500.0) & (wavelength < 1.0e5)
        kappa = np.loadtxt(dustdir+'dust01_trans.dat')
        w = w[mask]
        flux = flux[mask]
        kappa = kappa[mask]

    if sptype=='i' or sptype=='intrinsic':
        sp = S.ArraySpectrum(w, flux, fluxunits='flam')
        return(sp)

    a_dust = get_dust(w, p[0], model=dust_model)
    obsflux = flux * 10**(-0.4*a_dust)
    bb_scale = simps(flux-obsflux, w)

    bb = kappa * (w*1.0e-5)**-5 * 1.0/(np.exp(143843215.0/(w*p[3]))-1.0)
    # We also need to renormalize the bb flux.  The fraction of total luminosity
    # in the bb part of the spectrum is "scale".  So here we can simply multiply
    # by p[1] * FLUX_SCALE
    normalize_bb = simps(bb, w)
    bb = bb_scale * bb / normalize_bb

    return(simps(bb, w)/FLUX_SCALE)

# Get model, observed flux density from star
# Inputs are:
#      x = Name of filter function (provided in filters/ directory)
#      p = [p0 = tau_V from Kochanek et al. (2012),
#           p1 = Scaling constant for star flux (should be Lsol),
#           p2 = Temperature of input model star
#           p3 = temperature of dust]
def get_ext_bb(p, rsg_model='10', dust_model='g2', sptype='all', masked=True):

    w = wavelength
    flux = get_rsg(p[1], p[2], model=rsg_model)

    if masked and len(w)!=5157:
        mask = (wavelength > 3500.0) & (wavelength < 1.0e5)
        kappa = np.loadtxt(dustdir+'dust01_trans.dat')
        w = w[mask]
        flux = flux[mask]
        kappa = kappa[mask]

    if sptype=='i' or sptype=='intrinsic':
        sp = S.ArraySpectrum(w, flux, fluxunits='flam')
        return(sp)

    a_dust = get_dust(w, p[0], model=dust_model)
    obsflux = flux * 10**(-0.4*a_dust)
    bb_scale = simps(flux-obsflux, w)

    if sptype=='s' or sptype=='star':
        sp = S.ArraySpectrum(w, obsflux, fluxunits='flam')
        return(sp)

    bb = kappa * (w*1.0e-5)**-5 * 1.0/(np.exp(143843215.0/(w*p[3]))-1.0)
    # We also need to renormalize the bb flux.  The fraction of total luminosity
    # in the bb part of the spectrum is "scale".  So here we can simply multiply
    # by p[1] * FLUX_SCALE
    normalize_bb = simps(bb, w)
    bb = bb_scale * bb / normalize_bb

    if sptype=='b' or sptype=='bb' or sptype=='dust' or sptype=='blackbody':
        sp = S.ArraySpectrum(w, bb, fluxunits='flam')
    elif sptype=='ss' or sptype=='scaled_star':
        sp = S.ArraySpectrum(w, obsflux, fluxunits='flam')
    else:
        sp = S.ArraySpectrum(w, obsflux + bb, fluxunits='flam')
        mask = sp.flux < 0

    if sptype=='all':
        Lsol = simps(sp.flux, sp.wave)*4*np.pi*(10*3.08568025e18)**2/(3.839e33)
        try:
            assert np.abs((Lsol-p[1])/p[1]) < 3.0e-2
        except AssertionError:
            print(p)
            print(Lsol)
            print(np.abs((Lsol-p[1])/p[1]))
            print(rsg_model, dust_model, masked)
            sys.exit()

    return(sp)

def get_rv(p, rsg_model='10', dust_model='g2'):

    w = wavelength
    flux = get_rsg(p[1], p[2], model=rsg_model)

    mask = (wavelength > 3500.0) & (wavelength < 1.0e5)
    kappa = np.loadtxt(dustdir+'dust01_trans.dat')
    w = w[mask]
    flux = flux[mask]
    kappa = kappa[mask]

    a_dust = get_dust(w, p[0], model=dust_model)
    obsflux = flux * 10**(-0.4*a_dust)
    bb_scale = simps(flux-obsflux, w)

def get_new_rsg(p, rsg_model='2', dust_model='sil'):


    tau = p[0]
    lum = p[1]
    teff = p[2]
    tdust = p[3]

    file = f'data/interpolate/rsg_{dust_model}_s{rsg_model}_fullgrid.pkl'
    models = pickle.load(open(file, 'rb'))

    wavelengths = np.linspace(2000, 100000, 5000)
    flux = models((wavelengths, tau, lum, teff, tdust))

    sp = S.ArraySpectrum(wavelengths, flux, waveunits='angstrom', fluxunits='flam')

    return(sp)

