from astropy import units as u

def color(r,g,b):
    return((r/255.,g/255.,b/255., 1.0))

black=color(0,0,0)
red=color(255,0,0)
blue=color(10,0,255)
green=color(12,83,0)
magenta=color(204,0,204)
goldenrod=color(239,139,8)
orange=color(204,102,0)

pallette = [red,orange,goldenrod,green,blue,magenta]

SOLAR_LUM_CGS = u.solLum.to(u.erg/u.second)

# Convert normalized pysynphot.blackbody to solar luminosity flux units
BB_SCALE = 1.1483949e19

# Converts pickles spectrum normalized to 1 to solar luminosity flux units
PICKLES_CONST = 3.08397467e-7

# Constant to scale 4*pi*tauV/kappaV*R_in * R_out to Msun with kappaV in cm2/g,
# Rin and Rout in Rsun
DUST_BB_MASS = 2.4319771e-12

# Constant to scale 4*pi*tauV/kappV*R_in * v_wind to Msun/yr with kappaV in
# cm2/g, Rin in Rsun, vwind in km/s
DUST_BB_WIND = 1.10346114e-10

# 10km/s is default
RSG_V_WIND = 10.0

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
    'u':4.239,
    'g':3.303,
    'r':2.285,
    'i':1.698,
    'z':1.263,
    'PS1_GPC1_g':3.172,
    'PS1_GPC1_r':2.271,
    'PS1_GPC1_i':1.682,
    'PS1_GPC1_z':1.322,
    'PS1_GPC1_y':1.087,
    'WFPC2_F218W':7.760,
    'WFPC2_F300W':4.902,
    'WFPC2_F450W':3.410,
    'WFPC2_F555W':2.755,
    'WFPC2_F606W':2.415,
    'WFPC2_F702W':1.948,
    'WFPC2_F814W':1.549,
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
    'WFC3_UVIS_F555W':2.855,
    'WFC3_UVIS_F600LP':1.781,
    'WFC3_UVIS_F606W':2.488,
    'WFC3_UVIS_F625W':2.259,
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
    'ACS_WFC_F775W':1.629,
    'ACS_WFC_F814W':1.526,
    'ACS_WFC_F850LP':1.243,
    # These are calculated from Rieke & Lebofsky 1985
    'SPITZER_IRAC_CH1':0.058,
    'SPITZER_IRAC_CH2':0.023,
    'SPITZER_IRAC_CH3':0.021,
    'SPITZER_IRAC_CH4':0.020,
}
