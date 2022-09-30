"""
    SED Fitter by C.D. Kilpatrick

    Takes an input data table of photometry and compares to MIST stellar
    evolution tracks, blackbodies, Pickles stellar SEDs, BPASS models, and a
    red supergiant model from Kilpatrick & Foley (2018).
    Uses emcee to infer physical parameters of stellar systems from photometry.
"""
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

# Dependencies from this repo
import constants
import utilities
import dust

import pysynphot as S
S.setref(area = 25.0 * 10000)

# List of parameters to be fit for different types of SEDs
model_fit_params = {
    # Initial mass of primary, mass ratio, binary period
    'bpass': ['mass','ratio','period'],
    # Initial mass and stellar age
    'mist': ['mass','age'],
    'mist_terminal': ['mass'],
    'blackbody': ['luminosity','temperature'],
    'pickles': ['luminosity','temperature'],
    'rsg': ['tau_V', 'luminosity','temperature','dust_temp']
}

class sed_fitter(object):
    def __init__(self, verbose=True, interpolate=True):

        # Make sure dustmap is initialized
        utilities.import_dustmap()

        self.filename = ''

        self.usagestring = 'analysis.py objname model'

        # Handle file and directory structure for interpolated data and model
        # data from pickles, mist
        self.files = {
                'pickles': {'input': 'pickles.dat',
                            'interp': 'interpolate/pickles.pkl'},
                'rsg': {'interp': 'interpolate/rsg.pkl'},
                'mist': {'suffix': {'cmd': 'M.track.eep.cmd',
                                    'theoretical':'M.track.eep'}},
                'blackbody':{'interp':'interpolate/blackbody.pkl'},
                'extinction':{'file':'extinction'}
        }

        self.dirs = {
            'bandpass':'data/bandpass/',
            'extinction': 'data/extinction/',
            'input': 'data/input/',
            'mist': 'data/mist/',
            'bpass': 'data/bpass/',
            'figures': 'figures/',
            'backends': 'data/backends/',
            'data': 'data/'
        }

        for key in self.dirs.keys():
            if not os.path.exists(self.dirs[key]):
                os.makedirs(self.dirs[key])

        # Report output values to this many significant figures
        self.significant_figures = 4

        # For parameters that are prohibitive to calculate with every sample,
        # take this many random samples first and then calculate
        self.nsamples = 20000

        # Magnitude system of input data (AB is usually best)
        self.magsystem = 'abmag'

        # Storing model data and magnitudes for models
        self.models = None
        self.model_mags = {}

        self.model_size = 1000
        self.ages_shape = [0, 0, self.model_size]
        self.mass_shape = [0, 0, self.model_size]

        self.phottable = None
        self.backend = None
        self.pickles_data = None

        self.terminal = True    # Fit to terminal states of MIST and BPASS
        self.interpolate = interpolate # Fit to interpolated BB and Pickles SEDs
        self.limits = False     # Check during import if any values are limits
        self.extinction_model = False
        self.verbose = verbose

        self.options = None

        self.all_ages = []
        self.all_masses = []

        self.model_type = ''
        self.model_fit_params = []
        self.model_fit_blobs = []

        self.mist_masses = np.arange(7.5, 40.0, 0.5)
        self.metallicity = 0.014 # Metallicity in terms of Z

        # Distance and uncertainty in Mpc
        self.distance = [12.3, 1.8]
        self.dm = 5.0 * np.log10(self.distance[0]) + 25.0

        # For extinction likelihood values to use as a prior for modeling
        self.extinction = {
            'Av': None, 'Rv': None, 'pdf': None, 'cdf': None,
            'interpolation': None, 'likelihood': False,
            # For the extinction in a particular bandpass given Av,Rv
            'function': {}
        }

        # Coordinate for looking up MW extinction
        self.coord = None
        self.mw_ebv = 0.0

        # Can input host extinction values directly rather than use a model
        self.host_ext = None
        self.host_ext_inst_filt = {}

        # A list of rv values for input filters
        self.rv = []

        # Default wavelength binset in angstroms for pysynphot
        #self.waves = 10.0 + 10.0*np.arange(24000)
        self.waves = 3500.0 + 10.0*np.arange(9650)
        self.bandpasses = None

        self.file_suffix = {'theoretical':'M.track.eep',
                            'cmd': 'M.track.eep.cmd'}

        self.bounds = {
            'luminosity': [2.0, 7.0],
            'temperature': [2000., 100000.],
            'tau_V': [0.01, 6.0],
            'dust_temp': [200., 2000.],
            'mass': [0.1, 120.0],
            'age': [1.0e4, 13.0e9],
            'period': [0.0, 4.0], # in log10(days)
            'ratio': [0.1, 0.9],
            'Av': [0.0, 6.0],
            'Rv': [2.0, 6.0]
        }

        self.model_functions = {
            'blackbody': self.compute_blackbody_mag,
            'bpass': self.compute_bpass_mag,
            'pickles': self.compute_pickles_mag,
            'rsg': self.compute_rsg_mag,
        }

        self.load_model_functions = {
            'blackbody': self.load_blackbodies,
            'bpass': self.load_bpass,
            'mist': self.load_mist,
            'mist_terminal': self.load_mist,
            'pickles': self.load_pickles,
            'rsg': self.load_rsg,
        }

        if self.interpolate:
            self.model_functions['blackbody']=self.get_interpolated_mag
            self.model_functions['pickles']=self.get_interpolated_mag
            self.model_functions['rsg']=self.get_interpolated_mag

        # For BPASS or other models that use Vega mag units
        self.ab_to_vega = {}

    def add_options(self, parser=None, usage=None):
        import argparse
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,
                conflict_handler='resolve')
        parser.add_argument('objname', type=str,
            help='Object name to analyze.  Must have a phot file in input dir')
        parser.add_argument('model', type=str,
            help='Type of model to fit to object '+\
            '[blackbody|pickles|mist|bpass|rsg]')
        parser.add_argument('--niter','-n', type=int, default=1,
            help='Number of iterations over which to perform emcee')
        parser.add_argument('--nsteps', type=int, default=500,
            help='Number of steps to perform with emcee')
        parser.add_argument('--nwalkers', type=int, default=1000,
            help='Number of walkers to pass to emcee')
        parser.add_argument('--extinction', default=False, action='store_true',
            help='Use an input extinction file as opposed to assuming '+\
            'the host extinction value in the photometry file')
        parser.add_argument('--sigma', type=float, default=1.0,
            help='Vary the initial parameter guesses for each walker by'+\
            'best guess * np.random.lognormal(1.0, sigma)')
        parser.add_argument('--nsamples', type=int, default=20000,
            help='For rsg model, number of samples to use for calculating '+\
            'dust parameters.')
        parser.add_argument('--notau', default=False, action='store_true',
            help='For the RSG model, set this flag to ignore tau_V and Tdust '+\
            ' parameters as part of fit.')
        parser.add_argument('--skipdust', default=False, action='store_true',
            help='For the RSG model, set this flag to skip sampling the '+\
            'dust parameters as part of showing the output.')


        return(parser)

    def set_model_type(self, model_type, extinction=False, notau=False):
        self.model_type = model_type

        if model_type=='rsg' and notau:
            model_fit_params['rsg'] = ['luminosity','temperature']

        self.model_fit_params = model_fit_params[model_type]
        # Fit for host extinction with emcee blobs
        if extinction: self.model_fit_blobs = ['Av','Rv']

        # Restrict the temperture range if model type is RSG
        if model_type=='rsg':
            self.bounds['luminosity']=[3.0,5.6]
            self.bounds['temperature']=[2600.0, 8000.0]

        if model_type=='pickles':
            self.bounds['luminosity']=[3.0,6.0]
            self.bounds['temperature']=[2500.35,39810.7]

    # For a cumulative distribution function cdf with dimensions N, M
    # corresponding to arr1 and arr2, respectively, create a uniform random
    # variable and find the closest point in the cdf then return the
    # corresponding values from arr1 and arr2
    def inject_uniform_into_cdf(self, arr1, arr2, cdf, seed=None):
        if seed: val = seed
        else:
            val = np.random.uniform()
            val1 = np.random.uniform()
            val2 = np.random.uniform()
        sep1 = np.abs(arr1[1]-arr1[0]) ; sep2 = np.abs(arr2[1]-arr2[0])
        idx = np.unravel_index(np.argmin(np.abs(cdf-val)), cdf.shape)
        return(arr1[idx[1]]+2*(val1-0.5)*sep1, arr2[idx[0]]+2*(val2-0.5)*sep2)

    def load_extinction(self, fromfile='', val=None):

        extdir = self.dirs['extinction']
        # Check for host_ebv file name stored in phottable metadata
        if self.phottable and 'host_ebv' in self.phottable.meta.keys():
            if not utilities.is_number(self.phottable.meta['host_ebv']):
                fromfile = self.phottable.meta['host_ebv']
        # Files should be formatted following extinction ipynb
        if fromfile:
            self.extinction_model=True
            chifile = extdir+fromfile+'_chi.npy'
            avfile = extdir+fromfile+'_Av.npy'
            rvfile = extdir+fromfile+'_Rv.npy'
            pdf = np.load(chifile)
            pdf = pdf/np.min(pdf)
            cdf = 1./pdf

            av = np.load(avfile)
            rv = np.load(rvfile)

            self.extinction['cdf'] = cdf
            self.extinction['Av'] = av
            self.extinction['Rv'] = rv
            self.extinction['pdf'] = pdf

            self.extinction['interpolation'] = interpolate.interp2d(av, rv,
                pdf, kind='cubic', bounds_error=False)

            self.extinction['likelihood'] = True

        elif val:
            self.extinction_model=False
            self.host_ext=val # Input val must be formatted as (Av,Rv)
            mjd,inst_filt,mag,magerr = self.get_fit_parameters(self.phottable)
            bandpasses = [self.inst_filt_to_bandpass(i) for i in inst_filt]

            for i,bp in zip(inst_filt,bandpasses):
                self.host_ext_inst_filt[i]=self.calculate_extinction(val[0],
                    val[1], bp)

        elif (self.phottable and 'host_ebv' in self.phottable.meta.keys() and
            'host_rv' in self.phottable.meta.keys()):

            ebv = self.phottable.meta['host_ebv']
            rv = self.phottable.meta['host_rv']

            val = (ebv*rv, rv)

            self.extinction_model = False
            self.host_ext=val

            mjd,inst_filt,mag,magerr = self.get_fit_parameters(self.phottable)
            bandpasses = [self.inst_filt_to_bandpass(i) for i in inst_filt]

            for i,bp in zip(inst_filt,bandpasses):
                self.host_ext_inst_filt[i]=self.calculate_extinction(val[0],
                    val[1], bp)

        # Default just set host extinction to zero
        else:
            val = (0.0, 0.0)
            mjd,inst_filt,mag,magerr = self.get_fit_parameters(self.phottable)
            for i in inst_filt:
                self.host_ext_inst_filt[i]=0.0

    # Cardelli et al. extinction law.  Given input wavelength, Av, Rv, then
    # outputs a pysynphot spectrum with relative extinction versus wavelength
    def extinction_law(self, wave, Av, Rv, deredden=False):

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

    def get_ab_to_vega(self, bandpass):

        flat = np.zeros(len(self.waves))+1.0
        test = S.ArraySpectrum(self.waves, flat)
        kwargs = {'force': 'taper', 'binset': self.waves}
        obs = S.Observation(test, bandpass, **kwargs)

        # AB to Vega conversion is Mag_Vega - Mag_AB
        ab_to_vega = obs.effstim('vegamag')-obs.effstim('abmag')

        return(ab_to_vega)

    def import_phottable(self, filename):
        # Import phototmetry table - we require ascii.ecsv format
        try:
            table = Table.read(filename, format='ascii.ecsv')
        except astropy.io.ascii.core.InconsistentTableError:
            error = 'ERROR: input file needs to be ascii.ecsv format, e.g.,'
            print(error)
            sys.exit()

        # Mask the table to only include columns that we want to use
        if 'Use' in table.keys():
            mask = table['Use']==1
            table = table[mask]

        # Parse coord from table metadata
        if 'ra' in table.meta.keys() and 'dec' in table.meta.keys():
            self.coord = utilities.parse_coord(table.meta['ra'],
                table.meta['dec'])
            sfd = utilities.get_sfd()
            self.mw_ebv = sfd(self.coord)

        # Set mag system if it is in photometry table.  AB mag is preferred
        if 'magsystem' in table.meta.keys():
            system = table.meta['magsystem']
            if 'ab' in system.lower(): self.magsystem = 'abmag'
            if 'vega' in system.lower(): self.magsystem = 'vegamag'
            if 'st' in system.lower(): self.magsystem = 'stmag'

        if 'distance' in table.meta.keys():
            self.distance[0] = float(table.meta['distance'])
        if 'dist_error' in table.meta.keys():
            self.distance[1] = float(table.meta['dist_error'])

        self.phottable = table
        self.filename = filename

        # Rename phottable keys to lowercase
        for key in self.phottable.keys():
            self.phottable.rename_column(key, key.lower())

        # Rename common alternatives into standard column names
        for key in constants.alternate_names.keys():
            for val in constants.alternate_names[key]:
                if val in self.phottable.keys():
                    self.phottable.rename_column(val, key)

        # Check if we need to convert magnitudes to Vega mag
        if self.model_type=='bpass' and self.magsystem=='abmag':
            # Convert to Vega
            for i,row in enumerate(self.phottable):
                inst_filt = self.get_inst_filt(row)
                bp = self.inst_filt_to_bandpass(inst_filt)
                if bp.name not in self.ab_to_vega.keys():
                    self.ab_to_vega[bp.name] = self.get_ab_to_vega(bp)
                conv = self.ab_to_vega[bp.name]
                mag = self.phottable['magnitude'][i]
                self.phottable['magnitude'][i]=mag+conv

        # Check if we have any input values that are limits (error=0.0)
        error = list(self.phottable['mag_err'].data)
        if any([float(e)==0.0 for e in error]):
            self.limits = True

        # Set waveset depending on whether we're fitting Spitzer data
        wavemax = 2.5e4
        inst = list(self.phottable['instrument'].data)
        if any(['spitzer' in i.lower() for i in inst]): wavemax = 1.5e5
        S.setref(waveset=(100, wavemax, 50000, 'log'))

        # Sort table and get list of Rv values to apply for MW extinction
        self.phottable.sort('filter')
        newdict={}
        for key in constants.sf11rv.keys():
            newdict[key.upper()]=constants.sf11rv[key]
        constants.sf11rv.update(newdict)

        for row in self.phottable:
            inst_filt = self.get_inst_filt(row)
            self.rv.append(constants.sf11rv[inst_filt])

        return(self.phottable)

    # If we are using interpolated magnitudes, calculate an array of blackbody
    # magnitudes within the bounds on temperature and luminosity
    def load_blackbodies(self, phottable):

        inst_filt = [self.get_inst_filt(row) for row in phottable]
        # Make sure every value is unique
        inst_filt = list(np.unique(inst_filt))
        bandpasses = [self.inst_filt_to_bandpass(i) for i in inst_filt]

        # Luminosity is given in log-space, and calculate temperature with a
        # log distribution (although it is input as temp, not np.log10(temp))
        lum_bounds = self.bounds['luminosity']
        temp_bounds = np.log10(np.array(self.bounds['temperature']))
        Lval = np.linspace(lum_bounds[0], lum_bounds[1], 30)
        logTval = np.linspace(temp_bounds[0], temp_bounds[1], 50)

        if self.verbose:
            print('\n\nCreating grid of blackbodies with')
            print('logL_min={0}, logL_max={1}, T_min={2}, T_max={3}'.format(
                np.min(Lval), np.max(Lval), 10**np.min(logTval),
                10**np.max(logTval)))

        # Unpickle model mag function if it already exists
        bbfile = self.dirs['data']+self.files['blackbody']['interp']
        if os.path.exists(bbfile):
            models = pickle.load(open(bbfile, 'rb'))
            # Remove inst_filt pairs that we don't need to calculate
            for key in models.keys():
                if key in inst_filt:
                    if self.verbose:
                        print('Removing',key,'already in model file')
                    inst_filt.remove(key)
            # Exit if we have all models
            if not inst_filt:
                return(models)

        # Create models for the ones that don't already exist
        mags = {}
        for val in inst_filt:
            mags[val] = np.zeros((len(logTval), len(Lval)))

        bar = progressbar.ProgressBar(max_value=len(Lval)*len(logTval))
        bar.start()
        for i,L in enumerate(Lval):
            for j,logT in enumerate(logTval):
                mag = self.compute_blackbody_mag(inst_filt, L, 10**logT)
                bar.update(i*len(logTval)+j+1)
                for k,val in enumerate(inst_filt):
                    mags[val][j,i]=mag[k]

        bar.finish()

        models = {}
        for val in inst_filt:
            models[val] = interpolate.interp2d(Lval, logTval, mags[val],
                kind='cubic', bounds_error=True)

        # Save models back to pickle file
        pickle.dump(models, open(bbfile, 'wb'))

        return(models)

    # Load model magnitudes given the filters in the input phottable (i.e., for
    # passbands for photometry tables).  Magnitudes will be stored in an
    def load_models(self, phottable, **kwargs):

        if self.model_type not in self.load_model_functions.keys():
            error = 'ERROR: model type={type} is not allowed'
            print(error.format(type=self.model_type))

        func = self.load_model_functions[self.model_type]
        self.models = func(phottable, **kwargs)

    # Loads a file from the pickles stellar SED atlas.  Assumes pysynphot is
    # installed
    def import_pickles_file(self, filename):
        if '.fits' in filename: filename=filename.replace('.fits','')
        fullfile = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles',
            'dat_uvk', filename+'.fits')
        if os.path.exists(fullfile):
            sp = S.FileSpectrum(fullfile)
            return(sp)
        else:
            return(None)

    def load_pickles_data(self):

        # Import file with information about Pickles spectra
        filename = self.dirs['data']+self.files['pickles']['input']
        pickles_data = None
        if os.path.exists(filename):
            pickles_data = ascii.read(filename)
            for key in pickles_data.keys():
                pickles_data.rename_column(key, key.lower())
        if not pickles_data:
            error = 'ERROR: could not locate input pickles data'
            print(error)
            sys.exit()

        spectra = []
        for model in pickles_data:
            sp0 = self.import_pickles_file(model['filename'])
            sp0.convert('flam')

            # Normalize the spectrum so it integrates to 1 over range
            scale = 1./integrate.simps(sp0.flux, sp0.wave)
            sp0 = sp0 * scale * constants.PICKLES_CONST

            spectra.append(sp0)

        col = Column(spectra, name='spectrum')

        pickles_data.add_column(col)
        pickles_data.sort('temp')

        self.pickles_data = pickles_data

        return(pickles_data)

    def load_pickles(self, phottable, **kwargs):

        pickles_data = self.load_pickles_data()

        if not self.interpolate:
            return(pickles_data)
        else:
            # Calculate magnitudes for a range of luminosity values given
            # pickles spectra
            Lval = np.linspace(2.0, 7.0, 50)
            logTval = np.log10(pickles_data['temp'].data)
            LL, TT = np.meshgrid(Lval, logTval)
            LL = LL.ravel() ; TT = TT.ravel()

            inst_filt = [self.get_inst_filt(row) for row in phottable]
            remove_keys = []
            models = None

            pfile = self.dirs['data']+self.files['pickles']['interp']
            if os.path.exists(pfile):
                models = pickle.load(open(pfile, 'rb'))
                # Remove inst_filt pairs that we don't need to calculate
                for key in models.keys():
                    if key in inst_filt:
                        inst_filt.remove(key)
                # Exit if we have all models
                if not inst_filt:
                    return(models)

            # Create models for the ones that don't already exist
            mags = {}
            for val in inst_filt:
                mags[val] = np.zeros((len(logTval), len(Lval)))

            bandpasses = self.get_bandpasses(inst_filt)

            bar = progressbar.ProgressBar(max_value=len(Lval)*len(logTval))
            bar.start()
            for j,logT in enumerate(logTval):
                sp = pickles_data[j]['spectrum']
                for i,L in enumerate(Lval):
                    scaled_spec = sp * 10**L
                    mag = self.compute_pysynphot_mag(scaled_spec, bandpasses)
                    bar.update(j*len(Lval)+i+1)
                    for k,val in enumerate(inst_filt):
                        mags[val][j,i]=mag[k]

            bar.finish()

            models = {}
            for val in inst_filt:
                models[val] = interpolate.SmoothBivariateSpline(LL,
                    TT, list(mags[val].ravel()), kx=2,ky=2)

            # Save models back to pickle file
            pickle.dump(models, open(pfile, 'wb'))

            return(models)

    def load_rsg(self, phottable, **kwargs):

            # Calculate magnitudes for a range of luminosity, tau_V, and dust
            # temperatures
            #nLval = 4 ; nTdval = 4 ; ntau_Vval = 4 ; nTval = 4
            nLval = 26 ; nTdval = 12 ; ntau_Vval = 12 ; nTval = 20
            Lval = np.linspace(self.bounds['luminosity'][0],
                self.bounds['luminosity'][1], nLval)
            Tdval = 10**np.linspace(1, 3.35, nTdval)
            tau_Vval = 10**np.linspace(-2.0, 0.8, ntau_Vval)
            Tval = np.linspace(self.bounds['temperature'][0],
                self.bounds['temperature'][1], nTval)

            inst_filt = [self.get_inst_filt(row) for row in phottable]
            # Get unique inst_filt values
            inst_filt = list(np.unique(inst_filt))
            remove_keys = []
            models = None

            pfile = self.dirs['data']+self.files['rsg']['interp']
            models = {}
            if os.path.exists(pfile):
                models = pickle.load(open(pfile, 'rb'))
                # Remove inst_filt pairs that we don't need to calculate
                for key in models.keys():
                    if key in inst_filt:
                        inst_filt.remove(key)
                # Exit if we have all models
                if not inst_filt:
                    return(models)

            # Create models for the ones that don't already exist
            mags = {}
            for val in inst_filt:
                mags[val] = np.zeros((len(tau_Vval), len(Lval),
                    len(Tval), len(Tdval)))

            bandpasses = self.get_bandpasses(inst_filt)

            bar = progressbar.ProgressBar(max_value=len(Lval)*len(Tval)*\
                len(Tdval)*len(tau_Vval))
            bar.start()
            Nsteps=0
            for j, Teff in enumerate(Tval):
                for i,L in enumerate(Lval):
                    for k,Td in enumerate(Tdval):
                        for l,tau_V in enumerate(tau_Vval):
                            Nsteps=Nsteps+1
                            mag = self.compute_rsg_mag(inst_filt, tau_V, L,
                                Teff, Td)
                            bar.update(Nsteps)
                            for n,val in enumerate(inst_filt):
                                mags[val][l,i,j,k]=mag[n]

            bar.finish()

            params = (tau_Vval, 10**Lval, Tval, Tdval)
            for val in inst_filt:
                models[val] = interpolate.RegularGridInterpolator(params,
                    mags[val], method='linear', bounds_error=True)

            # Save models back to pickle file
            pickle.dump(models, open(pfile, 'wb'))

            return(models)

    def inst_filt_to_bandpass(self, inst_filt, bpdir='data/bandpass/'):
        inst = inst_filt.lower()
        filt = inst_filt.split('_')[-1].lower()

        bp = None
        if 'wfc3' in inst and 'uvis' in inst:
            if 'uvis2' in inst:
                bp=S.ObsBandpass('wfc3,uvis2,'+filt)
            # Default to UVIS1
            else:
                bp=S.ObsBandpass('wfc3,uvis1,'+filt)
        elif 'wfc3' in inst and 'ir' in inst:
            bp=S.ObsBandpass('wfc3,ir,'+filt)
        elif 'acs' in inst and 'wfc' in inst:
            if 'wfc2' in inst:
                bp=S.ObsBandpass('acs,wfc2,'+filt)
            # Default to WFC1
            else:
                bp=S.ObsBandpass('acs,wfc1,'+filt)
        elif 'acs' in inst and 'hrc' in inst:
            bp=S.ObsBandpass('acs,hrc,'+filt)
        elif 'acs' in inst and 'sbc' in inst:
            bp=S.ObsBandpass('acs,sbc,'+filt)
        elif 'wfpc2' in inst:
            bp=S.ObsBandpass('wfpc2,'+filt)
        else:
            # Load from bandpass directory
            if ('spitzer' in inst and 'irac' in inst):
                filename = 'SPITZER.IRAC.'+filt.upper()+'.dat'
            elif ('ps1' in inst and 'gpc1' in inst):
                filename = 'PS1.GPC1.'+filt.lower()+'.dat'
            elif ('lsst' in inst):
                filename = 'LSST.'+filt.lower()+'.dat'
            elif ('lcogt' in inst):
                filename = 'LCOGT.'+filt.lower()+'.dat'
            elif ('swift' in inst):
                filename = 'SWIFT.'+filt.lower()+'.dat'
            elif ('swope' in inst):
                filename = 'SWOPE.'+filt.lower()+'.dat'
            elif ('atlas' in inst):
                filename = 'ATLAS.'+filt.lower()+'.dat'
            elif ('asassn' in inst):
                filename = 'ASASSN.'+filt.lower()+'.dat'
            elif ('decam' in inst):
                filename = 'DECAM.'+filt.lower()+'.dat'

            file = bpdir + filename
            if not os.path.exists(file):
                print(f'WARNING: filter for {inst_filt} does not exist!')
                return(None)

            wave, trans = np.loadtxt(file, unpack=True)
            bp = S.ArrayBandpass(wave, trans, name=filt.lower(),
                waveunits='Angstrom')

        return(bp)

    # Given input list of inst_filt values, output a list of corresponding
    # bandpasses.  Output will be a list of lists, where elements in list
    # can correspond to a single bandpass or two for a color
    def get_bandpasses(self, inst_filt, redo=False, save=True):

        if (self.bandpasses and not redo and isinstance(inst_filt, list) and
            len(inst_filt)==len(self.bandpasses)):
            return(self.bandpasses)

        bpdir = self.dirs['bandpass']

        # Get bandpasses for each filter in phottable
        bandpasses = []
        for inst in inst_filt:
            if '-' in inst:
                vals = inst.split('-')
                bandpasses.append([self.inst_filt_to_bandpass(vals[0]),
                    self.inst_filt_to_bandpass(vals[1])])
            else:
                bandpasses.append([self.inst_filt_to_bandpass(inst)])

        if save:
            self.bandpasses = bandpasses

        return(bandpasses)

    def compute_pysynphot_mag(self, spectrum, bandpasses):
        # Get bandpasses for each filter in phottable
        mags = []
        kwargs = {'force': 'taper', 'binset': self.waves}
        for bp in bandpasses:

            if len(bp)==1:
                obs = S.Observation(spectrum, bp[0], **kwargs)

                try:
                    mag = obs.effstim(self.magsystem)
                    mags.append(mag)
                except ValueError:
                    mags.append(np.nan)

            elif len(bp)==2:
                obs1 = S.Observation(spectrum, bp[0], **kwargs)
                obs2 = S.Observation(spectrum, bp[1], **kwargs)

                try:
                    mag1 = obs1.effstim(self.magsystem)
                    mag2 = obs2.effstim(self.magsystem)
                    mags.append(mag1-mag2)
                except ValueError:
                    mags.append(np.nan)

        return(mags)

    def blackbody(self, temp, wave=None):
        if not wave:
            wave = 0.1 + 0.1*np.arange(9e5)

        flux = 1./wave**5 * 1./(np.exp(143843215./(wave*temp))-1.0)

        bb = S.spectrum.ArraySourceSpectrum(wave=wave, flux=flux,
            waveunits='angstrom', fluxunits='flam')

        return(bb)

    def create_blackbody(self, lum, temp):
        bb = S.BlackBody(temp)
        bb.convert('fnu')

        scale = constants.BB_SCALE/temp**4 * 10**lum
        bb = scale * bb

        return(bb)

    def compute_blackbody_mag(self, inst_filt, *args):

        # Get blackbody with right temperature and lum
        bb = self.create_blackbody(*args)

        # Get bandpasses and compute mags
        bandpasses = self.get_bandpasses(inst_filt)
        mags = self.compute_pysynphot_mag(bb, bandpasses)

        return(mags)

    def get_interpolated_mag(self, inst_filt, *args):
        if self.model_type=='pickles' or self.model_type=='blackbody':
            lum, temp = args
            mags = np.array([np.asscalar(self.models[i](lum, np.log10(temp)))
                for i in inst_filt])
        elif self.model_type=='rsg':
            if len(args)==2:
                L, Teff = args
                tau_V = 0.01
                Td = 200.
            else:
                tau_V, L, Teff, Td = args
            scaled_L = 10**L
            mags = np.array([np.asscalar(self.models[i]((tau_V, scaled_L,
                Teff, Td))) for i in inst_filt])
        return(mags)

    def create_pickles(self, lum, temp):
        if not self.pickles_data: self.pickles_data = self.load_pickles_data()
        idx = np.argmin(abs(self.pickles_data['temp']-temp))
        row = self.pickles_data[idx]

        # Scale spectrum up to input luminosity
        sp = row['spectrum'] * 10**(lum)
        sp.name = row['sptype']

        return(sp)

    def create_rsg(self, tau_V, lum, temp, dust_temp, sptype='all'):
        # RSG model takes luminosity in Lsol as input
        scaled_lum = 10**lum
        spec = dust.get_ext_bb((tau_V, scaled_lum, temp, dust_temp),
            sptype=sptype)

        return(spec)

    def compute_pickles_mag(self, inst_filt, lum, temp):

        # For temperature, pick the closest temperature in pickles data
        if ('temp' not in self.models.keys() or
            'spectrum' not in self.models.keys()):
            return(None)

        # Scale spectrum up to input luminosity
        sp = self.create_pickles(lum, temp)

        # Get bandpasses and compute mags
        bandpasses = self.get_bandpasses(inst_filt)
        mags = self.compute_pysynphot_mag(sp, bandpasses)

        return(mags)

    def compute_rsg_mag(self, inst_filt, tau_V, lum, temp, dust_temp):

        # Scale spectrum up to input luminosity
        sp = self.create_rsg(tau_V, lum, temp, dust_temp)

        # Get bandpasses and compute mags
        bandpasses = self.get_bandpasses(inst_filt)
        mags = self.compute_pysynphot_mag(sp, bandpasses)

        return(mags)

    def load_mist(self, phottable, mist_type='full', feh=0.0, rotation=False):

        # First parse photometry data into the format required by MIST table key
        photkeys = []
        for row in phottable:
            if inst=='WFC3/UVIS':
                inst_str = 'WFC3' ; inst_key = 'WFC3_UVIS'
            elif 'WFC3' in inst and 'IR' in inst:
                inst_str = 'WFC3' ; inst_key = 'WFC3_IR'
            elif 'ACS' in inst and 'WFC' in inst:
                inst_str = 'ACS_WFC' ; inst_key = 'ACS_WFC'
            elif 'ACS' in inst and 'HRC' in inst:
                inst_str = 'ACS_HRC' ; inst_key = 'ACS_HRC'
            elif 'WFPC2' in inst:
                inst_str = 'WFPC2' ; inst_key = 'WFPC2'
            else:
                warning = 'WARNING: could not interpret instrument string {0}'
                print(warning.format(inst))
                continue

            photkey = inst_key + '_' + filt.upper()
            photkeys.append((inst_str, photkey))

        # Get feh key
        feh_str = str(int(kwargs['feh']*100)).zfill(4)
        opt = self.files['mist']

        all_models = None

        for m in self.mist_masses:

            mass_str = str(int(m*10000)).zfill(7)
            model_table = None

            # So we only have to load each table once, iterate through the
            # unique instrument values
            for inst in np.unique([key[0] for key in photkeys]):

                directory = opt['dir'] + 'FEH_{feh}/{inst}/'
                directory = directory.format(feh=feh_str, inst=inst)
                fullfile = directory + mass_str + opt['suffix']['cmd']

                if not os.path.exists(fullfile):
                    warning = 'WARNING: data do not exist for:{0}'
                    print(warning.format(directory))
                    continue

                cmd_table = ascii.read(fullfile, header_start=14)

                if not model_table:
                    model_table = cmd_table['star_age','log_Teff','log_L']
                    masscol = Column([m]*len(model_table),name='mass')
                    fehcol = Column([feh]*len(model_table),name='feh')
                    model_table.add_column(masscol)
                    model_table.add_column(fehcol)

                # Load model data for each relevant photkey
                for photkey in photkeys:
                    if photkey[0]!=inst:
                        continue

                    if (photkey in cmd_table.keys() and
                       key not in model_table.keys()):

                        model_table.add_column(cmd_table[key])

                # Add full set of model data to self.models
                if not all_models:
                    all_models = model_table
                else:
                    all_models = vstack([all_models, model_table])

        # Now that we've assembled the models, interpolate with grid to obtain a
        # more precise estimate at a range of model parameters
        if mist_type=='terminal':
            ages = np.log10(self.models['star_age'].data)
            mass = self.models['mass'].data

            self.ages_shape = (np.min(ages), np.max(ages), self.model_size)
            self.mass_shape = (np.min(mass), np.max(mass), self.model_size)

            ages_vector = np.linspace(*self.ages_shape)
            mass_vector = np.linspace(*self.mass_shape)

            for row in phottable:
                inst = row['instrument'].upper().replace('/','_')
                filt = row['filter'].upper()
                inst_filt = inst+'_'+filt

                mass_input = []
                for m in mass:
                    # Interpolate in mass, but only take the terminal age
                    mass_val= self.models[self.models['mass']==m][-1][inst_filt]
                    mass_input.append(mass_val)

                mag_func = interpolate.interp1d(mass, mass_input)
                self.model_mags[inst_filt] = mag_func(mass_vector)

        else:
            ages = np.log10(self.models['star_age'].data)
            mass = self.models['mass'].data

            self.ages_shape = (np.min(ages), np.max(ages), self.model_size)
            self.mass_shape = (np.min(mass), np.max(mass), self.model_size)

            ages_vector = np.linspace(*self.ages_shape)
            mass_vector = np.linspace(*self.mass_shape)

            ages_grid, mass_grid = np.meshgrid(ages_vector, mass_vector)

            for row in phottable:
                inst = row['instrument'].upper().replace('/','_')
                filt = row['filter'].upper()
                inst_filt = inst+'_'+filt

                mag = interpolate.griddata((ages, mass),
                    self.models[inst_filt].data, (ages_grid, mass_grid),
                    method='cubic')

                self.model_mags[inst_filt] = mag

            # For all masses, if ages_vector > max(age for that mass), then set
            # the mag value to np.Infinity
            for i,m in enumerate(mass_vector):
                # Get closest mass in original table to this value
                close_mass = self.masses[np.argmin(abs(self.masses-m))]
                mass_table = self.models[self.models['mass']==close_mass]

                # Now get ages in ages_vector that are larger than maximum age
                # for mass_table
                max_age = np.log10(np.max(mass_table['star_age']))
                ages_idx = np.where(ages_vector > max_age)

                # If no idx to reset, skip
                if len(ages_idx[0])==0:
                    continue

                for row in phottable:
                    inst = row['instrument'].upper().replace('/','_')
                    filt = row['filter'].upper()
                    inst_filt = inst+'_'+filt

                    # For all model_mags, set to large value for ages_idx and
                    # mass value
                    self.model_mags[inst_filt][ages_idx,i] = None

        return(all_models)

    def load_bpass(self, phottable, bpass_type='terminal', feh=0.014):

        # First parse photometry data into the format required by MIST table key
        filts = list(phottable['filter'].data)

        # From p. 22 of BPASSv2.2.1_Manual.pdf
        # drive.google.com/drive/folders/1FvwVO3bBiswoMYsBueDrAhfMKs2By0rr
        filtmap = {'U':53,'B':54,'V':55,'R':56,'I':57,'J':58,'H':59,'K':60,
            'u':61,'g':62,'r':63,'i':64,'z':65,
            'f300w':66,'f336w':67,'f435w':68,'f450w':69,
            'f555w':70,'f606w':71,'f814w':72,
            'f625w':63,'f438w':68}

        feh_str = 'z'+str(int(feh*1000)).zfill(3)
        file = self.filename.split('.')[0]
        bpassfile = self.dirs['bpass']+'bpass_{0}_{1}.dat'.format(feh_str,
            bpass_type)

        table = None
        if os.path.exists(bpassfile):
            table = Table.read(bpassfile, format='ascii')

        photkeys = [] ; remove_filts = []
        for row in phottable:
            filt = row['filter'] ; inst = row['instrument']
            if inst=='WFC3/UVIS':
                inst_key = 'WFC3_UVIS'
            elif 'WFC3' in inst and 'IR' in inst:
                inst_key = 'WFC3_IR'
            elif 'ACS' in inst and 'WFC' in inst:
                inst_key = 'ACS_WFC'
            elif 'ACS' in inst and 'HRC' in inst:
                inst_key = 'ACS_HRC'
            elif 'WFPC2' in inst:
                inst_key = 'WFPC2'
            else:
                continue

            key = inst_key+'_'+filt

            if table and key in table.keys():
                remove_filts.append(filt)
                continue

            photkeys.append(key)

        for filt in remove_filts:
            filts.remove(filt)

        if not photkeys:
            return(table)

        elif table:
            for key in photkeys:
                newcol = Column([None]*len(table), name=key)
                table.add_column(newcol)

        # Get feh key
        pattern = self.dirs['bpass']+feh_str+'/sneplot-*'
        model_files = glob.glob(pattern)

        bar = progressbar.ProgressBar(max_value=len(model_files))
        bar.start()
        for i,model in enumerate(model_files):
            bar.update(i+1)

            cmd_table = Table.read(model, format='ascii')
            path, base = os.path.split(model)
            data = base.split('-')

            mass = float(data[2])
            ratio = float(data[3])
            period = float(data[4])

            add_row = True
            if table:
                mask = (table['mass']==mass) & (table['ratio']==ratio) &\
                    (table['period']==period)

                if len(table[mask])==1:
                    add_row = False
                    idx = np.where(mask)[0]
                    if bpass_type=='terminal':
                        row = cmd_table[-1]
                        for key,filt in zip(photkeys,filts):
                            if filt.lower() not in filtmap.keys():
                                continue
                            colnum = filtmap[filt.lower()]
                            mag = row['col'+str(colnum+1)]
                            table[idx][key]=mag

            if add_row:
                if bpass_type=='terminal':
                    row = cmd_table[-1]
                    table_data = [mass, ratio, period]

                    for filt in filts:
                        colnum = filtmap[filt.lower()]
                        mag = row['col'+str(colnum)]
                        table_data.append(mag)

                    if not table:
                        table_header = ['mass','ratio','period']+photkeys

                        table_data = [[dat] for dat in table_data]
                        table = Table(table_data, names=table_header)
                    else:
                        table.add_row(table_data)

        bar.finish()

        table.write(bpassfile, format='ascii', overwrite=True)

        return(table)

    def compute_bpass_mag(self, inst_filt, *args):
        # Split args into parameters for BPASS
        mass, ratio, period = args

        # Bisect mass then ratio then period to get closest row
        table = self.models

        for val in [('mass',mass),('ratio',ratio),('period',period)]:
            minimum = np.min(np.abs(table[val[0]]-val[1]))
            mask = np.abs(table[val[0]]-val[1])==minimum
            table = table[mask]

        # Should now be left with one row
        row = table[0]

        mags = []
        for data in inst_filt:
            if '-' in data:
                filt1,filt2 = data.split('-')
                mag1 = float(row[filt1]) ; mag2 = float(row[filt2])
                if mag1==99.0 or mag2==99.0:
                    mags.append(float('NaN'))
                else:
                    mags.append(mag1-mag2)
            elif data not in row.colnames:
                mags.append(float('NaN'))
            elif float(row[data])==99.0:
                mags.append(float('NaN'))
            else:
                mags.append(float(row[data]))

        return(mags)

    # Given input phottable row, output inst_filt variable
    def get_inst_filt(self, row):
        inst = row['instrument'].upper().replace('/','_')
        filt = row['filter'].upper()
        inst_filt = inst+'_'+filt
        return(inst_filt)

    # For a given phottable, outputs inst_filt, mag, magerr values that the
    # mcmc uses for fitting (as x, y, yerr, respectively)
    def get_fit_parameters(self, phottable):

        mjd = phottable['mjd'].data
        mag = phottable['magnitude'].data
        magerr = phottable['mag_err'].data
        inst_filt = np.array([r['instrument'].replace('/','_').upper()+'_'+\
            r['filter'].upper() for r in phottable])

        return(mjd, inst_filt, mag, magerr)

    def get_guess(self, model_type, guess_type='params'):

        if self.backend:
            try:
                if self.backend.iteration>0:
                    flat_samples = np.array(self.backend.get_chain(flat=True))
                    flat_prob = np.array(self.backend.get_log_prob(flat=True))
                    flat_blobs = np.array(self.backend.get_blobs(flat=True))

                    mask = np.isinf(np.abs(flat_prob))
                    flat_samples = flat_samples[~mask]
                    flat_prob = -1.0 * flat_prob[~mask]
                    flat_blobs = flat_blobs[~mask]

                    if len(flat_prob)>0:
                        best = np.argmin(flat_prob)
                        if guess_type=='params': return(flat_samples[best])
                        if guess_type=='blobs': return(flat_blobs[best])
            except OSError:
                pass

        if guess_type=='blobs': return(np.array([2.3, 4.4167]))

        guess = None
        if self.model_type=='mist_terminal': guess = np.array([30.0])
        elif self.model_type=='mist': guess = np.array([1.0e7, 30.0])
        elif self.model_type=='blackbody': guess = np.array([4.4, 3500.])
        elif self.model_type=='pickles': guess = np.array([5.077, 6353.])
        elif self.model_type=='bpass': guess = np.array([19.0, 0.1, 0.6])
        elif self.model_type=='rsg':
            if self.options.notau:
                guess = np.array([3.8, 3200])
            else:
                guess = np.array([0.02, 3.8, 3200, 1000])
        else:
            error = 'ERROR: unrecognized model type.  Exiting...'
            print(error)
            sys.exit()

        return(guess)

    def get_init_pos(self, guess, ndim, nwalkers, sigma=1, method='random'):

        init_pos = [guess * np.random.lognormal(1.0, sigma, ndim)
            for i in range(nwalkers)]

        return(init_pos)

    def load_backend(self, model_type, phottable, notau=False):
        # Formart objname from input filename
        objname = ''
        if 'name' in phottable.meta.keys():
            objname = phottable.meta['name']
        else:
            objname = self.filename
            if '/' in objname: objname = os.path.split(objname)[1]
            objname = objname.replace('.txt','')
            objname = objname.replace('.dat','')
            objname = objname.replace('.cat','')
            objname = objname.split('_')[0]
            objname = objname.split('-')[0]

        # Create a name from inst_filt, mag, magerr for comparison to backend
        mjd, inst_filt, mag, magerr = self.get_fit_parameters(phottable)
        name = ''
        for i,m,e in zip(inst_filt,mag,magerr):
            name += i+'='+str('%7.4f'%m)+'+/-'+str('%7.4f'%e)
            # Remove spaces
            name = name.replace(' ','')

        # I realized backend truncates names at 24 chars, so we need to generate
        # a (relatively) unique name limited to 24 chars.  I'm doing this by
        # taking every integer char in the name and modulo a 24 digit number.
        # This should result in the same name every time assuming the input
        # data (i.e., input bands, magnitudes, and errors) are the same.
        newname = ''
        for c in name:
            if c in '1234567890': newname+=c
        newname = str(int(newname)%100207100213100237100267)

        if self.model_type=='rsg' and notau:
            backfile = self.dirs['backends']+objname+'_rsg_notau.h5'
        else:
            backfile = self.dirs['backends']+objname+'_'+self.model_type+'.h5'
        if self.verbose:
            print('Backend file:',backfile)
            print('Backend name:',newname)
        backend = emcee.backends.HDFBackend(backfile, name=newname)

        return(backend)

    def sample_params(self, params, prob, ndim, nsamples=None, downsample=1.0):
        mask = np.isinf(np.abs(prob))
        if all(mask):
            print('WARNING: all probabilities are bad.  Try wider param range')
            return(params[0])
        if len(params.shape)==1:
            params = params[~mask]
        else:
            params = params[~mask,:]
        prob = -1.0 * prob[~mask]
        prob = prob / np.min(prob)

        # 1-sigma (0.6827 CL) is min chi^2 + 1.00 (1 param), 2.30 (2 param),
        # 3.50 (3 param), 4.72 (4 param), 5.89 (5 param), 7.04 (6 param)
        # https://ned.ipac.caltech.edu/level5/Wall2/Wal3_4.html provides the
        # values for 1-3 param.  The rest are calculated using Wolfram Alpha
        chi_limit = [1.00, 2.30, 3.50, 4.72, 5.89, 7.04]
        mask = prob < 1.0 + downsample * chi_limit[ndim-1]
        if len(params.shape)==1:
            params_sample = params[mask]
        else:
            params_sample = params[mask,:]
        prob_sample = prob[mask]

        if nsamples and nsamples < len(prob_sample):
            rand = np.array(random.sample(range(0, len(prob_sample)), nsamples))
            if len(params_sample.shape)==1:
                params_sample = params_sample[rand]
            else:
                params_sample = params_sample[rand,:]
            prob_sample = prob_sample[rand]

        return(params_sample, prob_sample)

    def calculate_param_best_fit(self, params, prob, ndim, name, show=True,
        sampled=False):

        # Parameters might have already been sampled
        if not sampled:
            params_sample, prob_sample = self.sample_params(params, prob, ndim)
        else:
            params_sample = params
            prob_sample = prob

        n = int(self.significant_figures)
        out_fmt = '{0:<18}: {1:>12} + {2:>12} - {3:>12}'

        # Get 16-50-84 values for param_sample now that we've restricted sample
        # to chi2 bounds
        best = np.percentile(params_sample, 50)
        minval = np.percentile(params_sample, 16)
        maxval = np.percentile(params_sample, 84)

        mcmc = utilities.round_to_n(best, n)
        digits = int(np.ceil(np.log10(mcmc)))
        decimal_place = -1 * (digits - n)
        if float(mcmc)==int(mcmc) and decimal_place < 1:
            mcmc=int(mcmc)

        minval = round(minval, decimal_place)
        maxval = round(maxval, decimal_place)
        maxval = maxval-mcmc
        minval = mcmc-minval

        if float(minval)==int(minval) and decimal_place < 1:
            minval=int(minval)
        if float(maxval)==int(maxval) and decimal_place < 1:
            maxval=int(maxval)

        if np.log10(mcmc)<-3:
            str_fmt = '%.3e'
            mcmc = str_fmt % mcmc
            maxval = str_fmt % maxval
            minval = str_fmt % minval
        elif decimal_place>0:
            str_fmt = '%7.{0}f'.format(int(decimal_place))
            mcmc = str_fmt % mcmc
            maxval = str_fmt % maxval
            minval = str_fmt % minval

        if show: print(out_fmt.format(name, mcmc, maxval, minval))

        return(best)

    def run_emcee(self, phottable, sigma=1.0, nsteps=5000, nwalkers=100,
        guess_type='params'):

        model_type = self.model_type

        # Indicate start of mcmc iteration
        self.make_banner(f'Starting MCMC for {model_type}')

        # This returns an absolute magnitude and colors for comparison for
        # the model data.  In this way, we can avoid over-estimating the
        # uncertainties due to distance and extinction by doing a filter-by-
        # filter comparison the absolute magnitudes
        mjd, inst_filt, mag, magerr = self.get_fit_parameters(phottable)

        # Get emcee parameters
        guess = self.get_guess(self.model_type, guess_type=guess_type)
        ndim = len(self.model_fit_params)

        # Load emcee backend
        if self.options: notau = self.options.notau
        backend = self.load_backend(self.model_type, self.phottable,
            notau=notau)

        # Decide if we're going to use an initial position from backend
        use_backend_pos = False
        if os.path.exists(backend.filename):
            try:
                if backend.iteration>0:
                    use_backend_pos = True
            except KeyError:
                pass

        # Check that nwalkers is consistent with the number of walkers on the
        # backend
        if use_backend_pos:
            walkers_backend = backend.shape[0]
            if nwalkers!=walkers_backend:
                m = 'WARNING: number of walkers on backend {0} does not '+\
                    'match input number of walkers {1}'
                m = m.format(walkers_backend, nwalkers)
                print(m)
                m = 'Do you want to change nwalkers to {0} [y/n]? '.format(
                    walkers_backend)
                resp = input(m)
                if resp=='y' or resp=='yes':
                    nwalkers = walkers_backend
                    self.options.nwalkers = walkers_backend
                else:
                    m = 'ERROR: inconsistent nwalkers on backend.\n'
                    m+= 'To solve this issue, use --nwalkers {0} '.format(
                        walkers_backend)
                    m+= 'or delete current backend.\n'
                    m+= 'Exiting...'
                    print(m)
                    sys.exit(1)

        if use_backend_pos:
            if self.verbose:
                print('Current number of iterations on backend:',
                    backend.iteration)
            init_pos = backend.get_last_sample()
        else:
            init_pos = self.get_init_pos(guess, ndim, nwalkers, sigma=sigma)

        # Construct emcee sampler with parameters derived above
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood,
            args=(inst_filt, mag, magerr), backend=backend,
            blobs_dtype=(float, float))

        # Run MCMC step - slow
        sampler.run_mcmc(init_pos, nsteps, progress=True)

        # Read current model probabilities, samples, and blobs from backend
        reader = self.load_backend(self.model_type, self.phottable)
        sample = np.array(reader.get_chain(flat=True))
        prob = np.array(reader.get_log_prob(flat=True))
        blob = np.array(reader.get_blobs(flat=True))

        if self.verbose:
            print('Current number of samples in backend:', len(sample))
            mask = ~np.isinf(prob)
            print('Minimum chi^2 is:','%.4f'%(-1.0*np.max(prob)))
            print('\n\n')
        params = [] ; blobs = []
        for i,param in enumerate(self.model_fit_params):
            p=self.calculate_param_best_fit(sample[:,i], prob, ndim, param)
            params.append(p)
        for i,param in enumerate(self.model_fit_blobs):
            b=self.calculate_param_best_fit(blob[:,i], prob, ndim, param)
            blobs.append(b)

        if self.distance[1]!=0.0 and self.verbose:
            print('Distance: ',self.distance[0],'+/-',self.distance[1])
            logL_unc = 2.17 * self.distance[1]/self.distance[0]
            print('Additional uncertainty on log_L:','%.4f'%logL_unc)
            print('Fractional uncertainty on log_R:','%.4f'%(logL_unc/2))

        # Calculate photospheric radius from best-fitting model params
        if ('luminosity' in self.model_fit_params and
            'temperature' in self.model_fit_params):

            lidx = np.array(['lum' in par for par in self.model_fit_params])
            tidx = np.array(['temperature' in par for par in self.model_fit_params])

            lum = 10**sample[:,lidx]
            temp = sample[:,tidx] / 5777
            radius = np.sqrt(lum / temp**4)

            r=self.calculate_param_best_fit(radius, prob, ndim, 'radius')

        # Calculate dust parameters from model params
        if (model_type=='rsg' and not self.options.skipdust and self.verbose
            and not self.options.notau):
            print('\n\n')
            print('RSG dust parameters:')
            # Need to sample parameters first since this part takes a long time
            # for a chain of samples > tens of thousands
            best_prob = np.array([np.max(prob)])
            extparams = np.array([params])
            param_sample, prob_sample = self.sample_params(sample, prob, ndim,
                nsamples=self.nsamples)
            param_sample = np.concatenate((param_sample, extparams), axis=0)
            prob_sample = np.concatenate((prob_sample, best_prob), axis=0)

            bar = progressbar.ProgressBar(max_value=len(param_sample))
            bar.start()
            bb_lum = []
            for i,p in enumerate(param_sample):
                # Account for the fact that log_L->L
                p[1] = 10**p[1]
                bar.update(i)
                bb_lum.append(dust.get_bb_lum(p))
            bar.finish()
            bb_lum = np.array(bb_lum)

            # First impose chi-2 cut and take random samples
            tidx = np.array(['dust_temp' in par
                for par in self.model_fit_params])
            kidx = np.array(['tau' in par for par in self.model_fit_params])

            dust_temp = param_sample[:,tidx] / 5777
            dust_temp = dust_temp[:,0]
            radius = np.sqrt(bb_lum / dust_temp**4)
            tau = param_sample[:,kidx]
            tau = tau[:,0]

            mass = 4 * np.pi * tau * dust.kappa_V**-1 * radius * 2 * radius *\
                constants.DUST_BB_MASS

            # Mass-loss rate assuming v_wind = 10 km/s
            mlr = 4 * np.pi * tau * dust.kappa_V**-1 * radius *\
                constants.DUST_BB_WIND * constants.RSG_V_WIND

            r=self.calculate_param_best_fit(np.log10(bb_lum), prob_sample, ndim,
                'dust luminosity', sampled=True)
            r=self.calculate_param_best_fit(radius, prob_sample, ndim,
                'dust radius', sampled=True)
            r=self.calculate_param_best_fit(mass, prob_sample, ndim,
                'dust mass', sampled=True)
            r=self.calculate_param_best_fit(mlr, prob_sample, ndim,
                'mass-loss rate', sampled=True)

        if self.verbose: print('\n\n')

        if not blobs and self.host_ext:
            blobs = self.host_ext

        model_mag, Av, Rv = self.compute_model_mag(inst_filt, params,
            extinction=blobs)

        if self.verbose:
            self.print_model_results(inst_filt, model_mag, mag, magerr, params,
                blobs, mjd)

    # Check bounds for input model parameters
    def check_bounds(self, theta):
        for i,par in enumerate(self.model_fit_params):
            if theta[i] < self.bounds[par][0] or theta[i] > self.bounds[par][1]:
                return(True)
        return(False)

    # Given an input set of inst_filts and model params theta, forward model to
    # the observed apparent mag (for given self.model_type, self.dm, extinction)
    def compute_model_mag(self, inst_filt, theta, extinction=None):

        model_mag = self.model_functions[self.model_type](inst_filt, *theta)

        # Catch bad model_mag value
        if model_mag is None:
            return(None, None, None)
        elif all([math.isnan(val) for val in model_mag]):
            return(None, None, None)

        # Apply dm and extinction according to probability distribution
        model_mag = np.array(model_mag) + self.dm
        if extinction is not None:
            Av, Rv = extinction
        else:
            Av, Rv = self.inject_uniform_into_cdf(self.extinction['Av'],
                self.extinction['Rv'], self.extinction['cdf'])

        for i,val in enumerate(inst_filt):
            if self.extinction_model:
                model_mag[i] += self.extinction['function'][val](Rv, Av)
            elif self.host_ext:
                model_mag[i] += self.host_ext_inst_filt[val]

            model_mag[i] += self.rv[i] * self.mw_ebv

        # Output model magnitudes and extinction values as blobs
        return(model_mag, Av, Rv)

    # Estimate log likelihood for a given age, mass, and data set
    def log_likelihood(self, theta, inst_filt, mag, magerr, extinction=None):

        if self.check_bounds(theta):
            return(-np.inf, None, None)

        if not extinction and self.host_ext:
            extinction=self.host_ext

        model_mag, Av, Rv = self.compute_model_mag(inst_filt, theta,
            extinction=extinction)

        # Catch bad model_mag value
        if model_mag is None:
            return(-np.inf, Av, Rv)

        # Flag missing data values
        mask = np.array(~np.isnan(model_mag))
        mag = mag[mask]
        magerr = magerr[mask]
        model_mag = model_mag[mask]

        # Handle limits
        if self.limits:
            limmask = magerr==0.0
            for m, mm in zip(mag[limmask], model_mag[limmask]):
                if m > mm:
                    return(-np.inf, Av, Rv)

            mag = mag[~limmask]
            magerr = magerr[~limmask]
            model_mag = model_mag[~limmask]

            # If all of the limits have passed check and there are no data
            # left, then we are in limit mode, so just return 1.0

            if len(mag)==0:
                return(-1.0, Av, Rv)

        chi2 = 1.0
        if self.extinction['likelihood']:
            chi2 *= self.extinction['interpolation'](Av, Rv)

        chi2 *= np.sum((mag-model_mag)**2/magerr**2)

        return(-1.0*chi2, Av, Rv)

    # Calculates extinction in input bandpass given Av and Rv
    def calculate_extinction(self, Av, Rv, bandpass):

        test = self.extinction_law(self.waves, Av, Rv)
        flat = np.zeros(len(self.waves))+1.0

        test1 = S.ArraySpectrum(self.waves, flat)
        test2 = S.ArraySpectrum(self.waves, flat*test.flux)

        kwargs = {'force': 'taper', 'binset': self.waves}

        obs1 = S.Observation(test1, bandpass, **kwargs)
        obs2 = S.Observation(test2, bandpass, **kwargs)

        # This is extinction in bandpass
        a_lambda = obs2.effstim('abmag')-obs1.effstim('abmag')

        return(a_lambda)

    # For a given set of bandpasses, calculate a grid of extinction values
    # A_filter given a range of extinction values Av and a range of reddening
    # values Rv
    def calculate_extinction_grid(self, inst_filt, Rv=np.linspace(1.0,6.0,30),
        Av=np.linspace(0.5,5.5,60)):

        bandpasses = [self.inst_filt_to_bandpass(i) for i in inst_filt]

        if self.verbose:
            print('\n\nCalculating extinction grid for',','.join(inst_filt))
            print('Using Av_min={0}, Av_max={1}, Rv_min={2}, Rv_max={3}'.format(
                np.min(Av), np.max(Av), np.min(Rv), np.max(Rv)))

        extinctionfile = self.dirs['extinction']+self.files['extinction']['file']
        if os.path.exists(extinctionfile+'.npy'):
            extinction_grid = np.load(extinctionfile+'.npy',
                allow_pickle=True)[()]
        else:
            extinction_grid = {}

        # Check that shape of grid is same as Av, Rv and that keys exist
        for val in inst_filt:
            if (val not in extinction_grid.keys() or
                extinction_grid[val].shape!=(len(Av), len(Rv))):
                extinction_grid[val] = None

        for i in np.arange(len(bandpasses)):
            bp = bandpasses[i] ; val = inst_filt[i]
            if val in extinction_grid.keys():
                if extinction_grid[val] is not None:
                    continue

            extinction_grid[val]=np.zeros((len(Av), len(Rv)))

            if self.verbose: print('Calculating extinction values for',val)

            for j,val1 in enumerate(Av):
                for k,val2 in enumerate(Rv):
                    outval = self.calculate_extinction(val1, val2, bp)
                    extinction_grid[val][j,k] = outval

        np.save(extinctionfile, extinction_grid, allow_pickle=True)

        # Interpolate the extinction grid so we can inject Av,Rv pairs and get
        # back an extinction value
        for key in extinction_grid.keys():
            if self.verbose: print('Running interpolation for',key)
            self.extinction['function'][key] = interpolate.interp2d(Rv, Av,
                extinction_grid[key], kind='cubic', bounds_error=False)

    def print_model_results(self, inst_filt, model, mag, magerr, theta, blobs,
         mjd, outtype=''):

        chi2, Av, Rv = self.log_likelihood(theta, inst_filt, mag, magerr,
            extinction=blobs)

        model_type = self.model_type
        if outtype=='initial':
            self.make_banner(f'Current best fit for: {model_type}')

        print('Input model parameters:')
        for name,par in zip(self.model_fit_params, theta):
            par = '%4.4f'%float(par)
            print(f'{name} = {par}')

        if self.model_fit_blobs:
            print('Input host extinction:')
            for name,par in zip(self.model_fit_blobs, blobs):
                par = '%4.5f'%par
                print(f'{name} = {par}')

        print('Av =', Av)
        print('Rv =', Rv)
        print('\n\n')

        # format output for each observation

        print('Model fit by observation:')
        fmt = '{mjd:<12} {inst:<20} {mag:<8} {err:<7} {model:<9} {chi:>10}'
        print(fmt.format(mjd='MJD', inst='INST_FILT',mag='MAG',err='MAGERR',
            model='MODEL',chi='CHI2'))
        for val in zip(inst_filt, mag, magerr, model, mjd):
            if float(val[2])!=0.0:
                chi=((val[1]-val[3])**2/val[2]**2)
                if chi>1.0e5:
                    chi='{:4e}'.format(chi)
                else:
                    chi='%5.4f'%chi
            elif val[3] < val[1]:
                chi='inf'
            else:
                chi='nan'
            print(fmt.format(mjd=val[4], inst=val[0], mag='%2.4f'%float(val[1]),
                err='%.4f'%float(val[2]), model='%2.4f'%float(val[3]),
                chi=chi))

        if self.extinction['likelihood']:
            value = np.asscalar(self.extinction['interpolation'](Av, Rv))
            print('Total Extinction chi^2:', value)

        if np.isinf(chi2):
            print('Sum chi^2 for best: inf')
        elif isinstance(chi2, float):
            print('Sum chi^2 for best: {0}'.format('%7.4f'%chi2))
        else:
            print('Sum chi^2 for best:', np.asscalar(-1.0*chi2))

    def get_wavelength_widths(self, spectrum, bandpasses):

        wave = [] ; width = []
        for bp in bandpasses:
            kwargs = {'force':'taper','binset': self.waves}
            obs = S.Observation(spectrum, bp, **kwargs)
            wave.append(obs.efflam())
            width.append(bp.rectwidth()/2.)

        return(np.array(wave), np.array(width))

    # Make sure all standard output is formatted in the same way with banner
    # messages for each module
    def make_banner(self, message):
        print('\n\n'+message+'\n'+'#'*80+'\n'+'#'*80+'\n\n')

    def parse_photfile(self, objname):
        filename = objname+'.dat'
        if os.path.exists(filename):
            return(filename)
        elif os.path.exists(self.dirs['input']+filename):
            return(self.dirs['input']+filename)
        else:
            return(None)

def main():
    # Start timer, create sed_fitter instance
    start = time.time()
    sed = sed_fitter()

    # Handle the --help option
    if '-h' in sys.argv or '--help' in sys.argv:
        parser = sed.add_options(usage=sed.usagestring)
        options = parser.parse_args()
        sys.exit()

    # Starting banner
    sed.command = ' '.join(sys.argv)
    sed.make_banner('Starting: {cmd}'.format(cmd=sed.command))

    parser = sed.add_options(usage=sed.usagestring)
    opt = parser.parse_args()
    sed.options = opt

    sed.nsamples = opt.nsamples

    photfile = sed.parse_photfile(opt.objname)
    if not photfile or not os.path.exists(photfile):
        inpfile = sys.argv[1]
        print(f'ERROR: input photfile {inpfile} does not exist!  Exiting...')
        sys.exit(1)

    if opt.model not in list(model_fit_params.keys()):
        model = opt.model
        print(f'ERROR: model type {model} not allowed!  Exiting...')
        sys.exit(1)

    run_models = [opt.model]*opt.niter
    for model in run_models:
        # Grab options again in case they were changed as part of an iteration
        opt = sed.options

        # Set model type and import photometry table, bandpasses
        # Input host extinction in (Av,Rv)
        sed.set_model_type(model, extinction=opt.extinction, notau=opt.notau)
        sed.phottable = sed.import_phottable(photfile)
        mjd, inst_filt, mag, magerr = sed.get_fit_parameters(sed.phottable)
        bandpasses = [sed.inst_filt_to_bandpass(i) for i in inst_filt]

        # Load extinction priors
        # This can be done with input E(B-V) and Rv values for the host or by
        # calculating a grid of extinction values for the input model parameters
        sed.load_extinction()
        #sed.calculate_extinction_grid(inst_filt)
        sed.load_models(sed.phottable)

        # Check if there is a backend for the current photometry table and model
        # and import or make one if there's not a current one
        sed.backend = sed.load_backend(model, sed.phottable, notau=opt.notau)

        # Get an initial guess for the model parameter values and compute the
        # model magnitudes and chi^2 for those model parameters
        theta = sed.get_guess(model, guess_type='params')
        ext = sed.host_ext
        model, Av, Rv = sed.compute_model_mag(inst_filt, theta, extinction=ext)
        chi2, Av, Rv = sed.log_likelihood(theta, inst_filt, mag, magerr,
            extinction=ext)

        if sed.verbose:
            sed.print_model_results(inst_filt, model, mag, magerr, theta, ext,
                mjd, outtype='initial')

        # Run the MCMC
        sed.run_emcee(sed.phottable, sigma=opt.sigma, nsteps=opt.nsteps,
            nwalkers=opt.nwalkers)

if __name__ == "__main__":
    main()
