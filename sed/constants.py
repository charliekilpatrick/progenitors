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
