import warnings
warnings.filterwarnings('ignore')

from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import vstack, Column, Table

from scipy import integrate, interpolate

import emcee
import os
import sys
import math
import glob
import pickle
import astropy
import progressbar
import numpy as np
import argparse
import time
import random

import pysynphot as S
S.setref(area = 25.0 * 10000)


waves = 100.0 + 10.0*np.arange(9650)

def extinction_law(wave, Av, Rv, deredden=False):

        # Inverse wavelength dependent quantities
        x=1./(wave*1.0e-4)
        y=x-1.82

        # Variables proportional to extinction
        a=np.zeros(len(wave))
        b=np.zeros(len(wave))
        Fa=np.zeros(len(wave))
        Fb=np.zeros(len(wave))

        # Masks to apply piecewise Cardelli et al. function
        mask1 = np.where(x < 1.1)
        mask2 = np.where((x > 1.1) & (x < 3.3))
        mask3 = np.where(x > 3.3)
        mask4 = np.where(x > 5.9)

        x=1./(wave[mask1]*1.0e-4)
        y=x-1.82

        a[mask1]=0.574 * x**1.61
        b[mask1]=-0.527 * x**1.61

        x=1./(wave[mask2]*1.0e-4)
        y=x-1.82

        a[mask2]=1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 +\
            0.72085 * y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b[mask2]=1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 -\
            5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7

        x=1./(wave[mask4]*1.0e-4)
        y=x-1.82

        Fa[mask4]=-0.04473 * (x-5.9)**2 - 0.009779 * (x-5.9)**3
        Fb[mask4]=0.2130 * (x-5.9)**2 + 0.1207 * (x-5.9)**3

        x=1./(wave[mask3]*1.0e-4)
        y=x-1.82

        a[mask3]=1.752 - 0.316 * x - 0.104/((x-4.67)**2 + 0.341) + Fa[mask3]
        b[mask3]=-3.090 + 1.825 * x + 1.206/((x-4.62)**2 + 0.263) + Fb[mask3]

        Alam = Av*(a+b/Rv)
        elam = 10**(-0.4 * Alam)
        if deredden: elam = 1./elam

        sp = S.ArraySpectrum(wave, elam, fluxunits='count')

        return(sp)

def calculate_extinction(Av, Rv, bandpass):

        test = extinction_law(waves, Av, Rv)
        flat = np.zeros(len(waves))+1.0

        test1 = S.ArraySpectrum(waves, flat)
        test2 = S.ArraySpectrum(waves, flat*test.flux)

        kwargs = {'force': 'taper', 'binset': waves}

        obs1 = S.Observation(test1, bandpass, **kwargs)
        obs2 = S.Observation(test2, bandpass, **kwargs)

        # This is extinction in bandpass
        a_lambda = obs2.effstim('abmag')-obs1.effstim('abmag')

        a_lambda = a_lambda * 0.88090941898

        return(a_lambda)


bp = S.ObsBandpass(sys.argv[1])
a_lambda = calculate_extinction(3.1, 3.1, bp)
print(a_lambda)
