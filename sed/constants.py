from astropy import units as u
from astropy.constants import astropyconst20 as const

def color(r,g,b):
    return((r/255.,g/255.,b/255., 1.0))

black=color(0,0,0)
red=color(255,0,0)
blue=color(10,0,255)
green=color(12,83,0)
magenta=color(204,0,204)
goldenrod=color(239,139,8)
orange=color(204,102,0)
lightred=color(255,178,178)

palette = [red,orange,goldenrod,green,blue,magenta]

SOLAR_LUM_CGS = u.solLum.to(u.erg/u.second)

# Convert normalized pysynphot.blackbody to solar luminosity flux units
BB_SCALE = 1.1483949e19

# Converts pickles spectrum normalized to 1 to solar luminosity flux units
PICKLES_CONST = 3.08397467e-7

# Constant to scale 4*pi*tauV/kappaV*R_in * R_out to Msun with kappaV in cm2/g,
# Rin and Rout in Rsun.  See run_emcee in analysis.py.
# This is essentially R_sun^2 * 1 cm^2/gram / (1 M_sun)
DUST_BB_MASS = 2.4319771e-12

# Constant to scale 4*pi*tau_V/kappa_V*R_in * v_wind to Msun/yr with kappa_V in
# cm2/g, Rin in Rsun, vwind in km/s.
# This is essentially (R_sun * 1 km/s) / (1 Msun/yr * 1 cm2/g)
d=const.R_sun.to('cm') * 1.0 * u.km/u.s / (const.M_sun.to('g') / (1.0 * u.year) * 1.0 * u.cm**2/u.g)
DUST_BB_WIND = d.to(u.Unit(1)).value

# 50km/s is default
RSG_V_WIND = 50.0 * u.km/u.s

FIGSIZE = 10.0

# Alternate names for phottable on input
alternate_names = {
    'magnitude': ['mag','m'],
    'mag_err': ['magerr','err','error','uncertainty'],
    'filter': ['band','bandpass','filt'],
    'instrument': ['inst','detector','det']
}

color_palette = {
    'WFC3_UVIS_F438W': blue,
    'WFC3_UVIS_F555W': green,
    'WFC3_UVIS_F625W': goldenrod,
    'WFC3_UVIS_F814W': red
}

# From Table 6 of Schlafly & Finkbeiner for R_filter values given Rv=3.1
sf11rv = {
    'U':4.334,
    'B':3.626,
    'V':2.742,
    'R':2.169,
    'I':1.505,
    'J':0.709,
    'H':0.449,
    'K':0.302,
    'KCONT':0.302,
    'KS':0.302,
    'u':4.239,
    'g':3.303,
    'r':2.285,
    'i':1.698,
    'z':1.263,
    'Y':1.2245,
    'ATLAS_O': 1.992,
    'ATLAS_C': 2.794,
    'PS1_GPC1_g':3.172,
    'PS1_GPC1_r':2.271,
    'PS1_GPC1_i':1.682,
    'PS1_GPC1_z':1.322,
    'PS1_GPC1_y':1.087,
    'DECAM_g': 3.237,
    'DECAM_Y': 1.058,
    'DECAM_z': 1.217,
    'DECAM_i': 1.595,
    'DECAM_r': 2.176,
    'WFPC2_F218W':7.760,
    'WFPC2_F300W':4.902,
    'WFPC2_F336W':4.453,
    'WFPC2_F450W':3.410,
    'WFPC2_F547M':2.801,
    'WFPC2_F555W':2.755,
    'WFPC2_F606W':2.415,
    'WFPC2_F656N':2.102,
    'WFPC2_F675W':2.169,
    'WFPC2_F702W':1.948,
    'WFPC2_F814W':1.549,
    'WFC3_IR_F098M':1.189,
    'WFC3_IR_F105W':0.969,
    'WFC3_IR_F110W':0.881,
    'WFC3_IR_F125W':0.726,
    'WFC3_IR_F140W':0.613,
    'WFC3_IR_F160W':0.512,
    'WFC3_UVIS_F200LP':2.958,
    'WFC3_UVIS_F218W':7.760,
    'WFC3_UVIS_F225W':6.989,
    'WFC3_UVIS_F275W':5.487,
    'WFC3_UVIS_F300X':5.228,
    'WFC3_UVIS_F336W':4.453,
    'WFC3_UVIS_F350LP':2.624,
    'WFC3_UVIS_F390W':3.896,
    'WFC3_UVIS_F438W':3.623,
    'WFC3_UVIS_F475W':3.248,
    'WFC3_UVIS_F475X':3.116,
    'WFC3_UVIS_F547M':2.801,
    'WFC3_UVIS_F555W':2.855,
    'WFC3_UVIS_F600LP':1.781,
    'WFC3_UVIS_F606W':2.488,
    'WFC3_UVIS_F621M':2.211,
    'WFC3_UVIS_F625W':2.259,
    'WFC3_UVIS_F673N':2.091,
    'WFC3_UVIS_F775W':1.643,
    'WFC3_UVIS_F814W':1.536,
    'WFC3_UVIS_F850LP':1.208,
    'ACS_CLEAR':2.436,
    'ACS_WFC_F435W':3.610,
    'ACS_WFC_F475W':3.268,
    'ACS_WFC_F550M':2.620,
    'ACS_WFC_F555W':2.792,
    'ACS_WFC_F606W':2.471,
    'ACS_WFC_F625W':2.219,
    'ACS_WFC_F658N':2.215,
    'ACS_WFC_F775W':1.629,
    'ACS_WFC_F814W':1.526,
    'ACS_WFC_F850LP':1.243,
    # These are calculated from Rieke & Lebofsky 1985
    'SPITZER_IRAC_CH1':0.058,
    'SPITZER_IRAC_CH2':0.023,
    'SPITZER_IRAC_CH3':0.021,
    'SPITZER_IRAC_CH4':0.020,
    # From https://iopscience.iop.org/article/10.3847/1538-4357/ab1c61/pdf
    'JWST_NIRCAM_F070W': 0.6919*3.1,
    'JWST_NIRCAM_F090W': 0.4523*3.1,
    'JWST_NIRCAM_F115W': 0.2785*3.1,
    'JWST_NIRCAM_F150W': 0.1618*3.1,
    'JWST_NIRCAM_F200W': 0.0919*3.1,
    'JWST_NIRCAM_F277W': 0.0470*3.1,
    'JWST_NIRCAM_F356W': 0.0283*3.1,
    'JWST_NIRCAM_F444W': 0.0194*3.1,
    'JWST_NIRCAM_F150W2': 0.1519*3.1,
    'JWST_NIRCAM_F322W2': 0.0370*3.1,
    # These come from https://arxiv.org/pdf/2303.04820.pdf
    'JWST_NIRCAM_F140M':0.315*3.1,
    'JWST_NIRCAM_F162M':0.255*3.1,
    'JWST_NIRCAM_F164N':0.251*3.1,
    'JWST_NIRCAM_F182M':0.214*3.1,
    'JWST_NIRCAM_F187N':0.209*3.1,
    'JWST_NIRCAM_F210M':0.180*3.1,
    'JWST_NIRCAM_F212N':0.176*3.1,
    'JWST_NIRCAM_F250M':0.143*3.1,
    'JWST_NIRCAM_F300M':0.116*3.1,
    'JWST_NIRCAM_F323N':0.107*3.1,
    'JWST_NIRCAM_F335M':0.103*3.1,
    'JWST_NIRCAM_F360M':0.096*3.1,
    'JWST_NIRCAM_F405N':0.088*3.1,
    'JWST_NIRCAM_F410M':0.087*3.1,
    'JWST_NIRCAM_F430M':0.084*3.1,
    'JWST_NIRCAM_F460M':0.080*3.1,
    'JWST_NIRCAM_F466N':0.079*3.1,
    'JWST_NIRCAM_F470N':0.079*3.1,
    'JWST_NIRCAM_F480M':0.078*3.1,
}
