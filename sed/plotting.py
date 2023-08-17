"""
    Methods for plotting SEDs and parameters inferred by sed_fitter
"""
import warnings
warnings.filterwarnings('ignore')

from astropy.io import ascii
from astropy.table import vstack, Column, Table
import os, numpy as np
from astropy.time import Time
import subprocess
from numpy.random import rand
import copy

from scipy import integrate

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from matplotlib.colors import ListedColormap

import corner

from analysis import *
from constants import *

import pysynphot as S
S.setref(area = 25.0 * 10000)

MPC_TO_CM = 3.08568025e24

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

class sed_plot(object):
    def __init__(self):

        self.pad = 0.5
        self.figsize = 12.0

        self.fit = sed_fitter(verbose=False)
        self.options = None

        self.outdir = 'figures/'

        rc('font',**{'family':'serif','serif':['Times'],
            'size':5.0*self.figsize})
        rc('text', usetex=True, color=black)

        flux_fmt = r'$_{\lambda}$'+' [10'+r'$^{0}$'
        flux = '{0}{1} {2}]'

        self.axis_titles = {
            'wavelength': 'Observer-frame Wavelength ('+r'\AA'+')',
            'time_mjd': 'Modified Julian Date',
            'time': 'Observer-frame Days from Discovery',
            'flux': 'L'+r'$_{\lambda}$'+' [10'+r'$^{23}$ {0}]',
            'L_lambda': flux.format('L',flux_fmt,
                ' erg s'+r'$^{-1}$~\AA$^{-1}$'),
            'f_lambda': flux.format('f', flux_fmt,
                ' erg s'+r'$^{-1}$'+' cm'+r'$^{-2}$~\AA$^{-1}$'),
            'log_T': r'$\log(T_{\mathrm{eff}}/\mathrm{K})$',
            'Teff': r'$T_{\mathrm{eff}}/\mathrm{K}$',
            'Tdust': r'$T_{\mathrm{dust,eff}}/\mathrm{K}$',
            'A_V': r'$A_{V}$',
            'log_L': r'$\log(L/L_{\odot})$',
            'initial_mass': r'$M/M_{\odot}$',
            'period': r'$\log(P/\mathrm{1 day})$',
            'mass_ratio': 'q',
        }

        self.plot_types = {
            'hr': self.hr,
            'sed': self.plot_sed,
            'cmd': self.cmd,
            'corner': self.plot_corner,
            'lightcurve': self.lightcurve,
        }

        self.usagestring = 'plotting.py objname type models'

    def add_options(self, parser=None, usage=None):
        import argparse
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage,
                conflict_handler='resolve')
        parser.add_argument('objname', type=str,
            help='Object name to analyze.  Must have a phot file in input dir')
        parser.add_argument('type', type=str,
            help='Types of models to show, can be comma-separated to show '+\
            'multiple models [blackbody|pickles|mist|bpass|rsg]')
        parser.add_argument('model', type=str,
            help='Types of plot to make, can be comma-separated to make '+\
            'multiple plots [hr|sed|cmd|corner|lightcurve]')
        parser.add_argument('--notau', default=False, action='store_true',
            help='For the RSG model, set this flag to ignore tau_V and Tdust '+\
            ' parameters as part of fit.')
        parser.add_argument('--simple-rsg', default=False, action='store_true',
            help='Do not add additional SEDs for the RSG model plot.')

        return(parser)

    def make_banner(self, message):
        print('\n\n'+message+'\n'+'#'*80+'\n'+'#'*80+'\n\n')

    def create_flux_title(self, power):
        flux_fmt = r'$_{\lambda}$'+' [10'+r'$^{'+str(power)+'}$'
        flux = '{0}{1} {2}]'

        title = flux.format('f', flux_fmt,
            ' erg s'+r'$^{-1}$'+' cm'+r'$^{-2}$~\AA$^{-1}$')
        return(title)

    def load_sed(self, model, photfile, file='', extinction=False,
        notau=False):

        # Set model type and import photometry table, bandpasses
        self.fit.set_model_type(model, extinction=extinction)
        self.fit.phottable = self.fit.import_phottable(photfile)
        mjd,inst_filt,mag,magerr = self.fit.get_fit_parameters(self.fit.phottable)
        bandpasses = [self.fit.inst_filt_to_bandpass(i) for i in inst_filt]

        # Load extinction priors
        self.fit.load_extinction(fromfile=file)

        # Calculate a grid of extinction values for the input model parameters
        if extinction: self.fit.calculate_extinction_grid(inst_filt)
        self.fit.load_models(self.fit.phottable)

        # Check if there is a backend for the current photometry table and model
        # and import or make one if there's not a current one
        self.fit.backend = self.fit.load_backend(model, self.fit.phottable,
            notau=notau)

    def convert_to_fluxspace(self, spectrum, ext):

        # Spectra should be scaled to absolute magnitudes, so flux at 10 pc,
        # so need to rescale by distance/10pc
        distance = self.fit.distance[0]
                
        # Assume distance in Mpc, so D/10pc = distance * 1e5
        distance = distance * 1e5
                
        # This is for host extinction
        extinction1 = self.fit.extinction_law(spectrum.wave, ext[0], ext[1])
                
        # This is for MW extinction
        Av = self.fit.mw_ebv * 3.1
        extinction2 = self.fit.extinction_law(spectrum.wave, Av, 3.1)
                
        # Now calculate observed flux and re-generate spectrum
        obsflux = spectrum.flux * extinction1.flux * extinction2.flux
        obsflux = obsflux/(distance**2)
        sp = S.ArraySpectrum(spectrum.wave, obsflux)

        return(sp)

    def plot_sed(self, models, fluxspace=True, extinction=False,
        notau=False, **kwargs):

        phot = self.fit.phottable
        dum, inst_filt, mag, magerr = self.fit.get_fit_parameters(phot)

        inst_filt = list(inst_filt)
        mag = list(mag)
        magerr = list(magerr)

        inst_filt = np.array(inst_filt)
        mag = np.array(mag)
        magerr=np.array(magerr)

        # Get flux for observations in Janskies
        flux = 3631 * 10**(-0.4 * np.array(mag)) * 1.0e-23
        fluxerr = flux * np.array(magerr)/1.086
        fluxcorr = False

        bandpasses = np.array([self.fit.inst_filt_to_bandpass(i) for i in inst_filt])

        # Set up plot
        fig, ax = self.setup_plot()
        self.setup_ticks(ax)

        for i,model in enumerate(models):

            print('Plotting model:',model)
            self.load_sed(model, self.fit.filename, extinction=extinction,
                notau=notau)

            self.fit.backend = self.fit.load_backend(model, self.fit.phottable,
                notau=notau)

            params = self.fit.get_guess(model, guess_type='params')
            ext = self.fit.get_guess(model, guess_type='blobs')

            sp = None
            if model=='blackbody':
                sp = self.fit.create_blackbody(*params)
                plottemp = ''
                title=str(int(np.round(params[1], decimals=-1)))+' K blackbody'
                sp.convert('flam')
            elif model=='pickles':
                self.fit.models = self.fit.load_pickles(self.fit.phottable)
                sp = self.fit.create_pickles(*params)
                title=sp.name
                title = title.replace('V','I')
                sp.convert('flam')
                sp = S.ArraySpectrum(sp.wave, sp.flux/1.086)
            elif model=='rsg':
                if len(params)==2:
                    title = 'RSG'
                    params = [0.01,params[0],params[1],200.]
                else:
                    label = 'Reddened RSG +\n {0} K dust'
                    dust_temp = int(truncate(params[3],-1))
                    title = label.format(dust_temp)

                sp = self.fit.create_rsg(*params)

            if fluxspace: sp = self.convert_to_fluxspace(sp, ext)

            # Get wavelengths for inst_filt and spectrum
            wave, width = self.fit.get_wavelength_widths(sp, bandpasses)
            # Convert flux (currently in erg/s/cm2/Hz) to flam
            # Only need to do this once
            if not fluxcorr:
                fluxnu = copy.copy(flux)
                flux = fluxnu * 2.998e18/wave**2
                fluxerr = fluxerr * 2.998e18/wave**2
                fluxcorr = True

            if i==0:
                print('filter mag wave flux_nu flux_lam')
                for j in np.arange(len(wave)):
                    if fluxerr[j]==0.0:
                        uplim = np.array([1])
                        fluxerr[j] = 0.3 * flux[j]
                        use_color=lightred
                    else:
                        uplim = np.array([0])
                        use_color=red

                    if magerr[j]!=0.0:
                        print(inst_filt[j],mag[j],wave[j],fluxnu[j],flux[j])
                    ax.errorbar([wave[j]], [flux[j]], yerr=[fluxerr[j]],
                        xerr=[width[j]], marker='o', color=use_color, ls='none',
                        capsize=0.5*self.figsize, linewidth=0.6*self.figsize,
                        zorder=5, ms=2*self.figsize, markeredgecolor=black,
                        markeredgewidth=0.4*self.figsize, uplims=uplim)

            if model=='blackbody':

                bbwave = sp.wave 
                bbsflux = sp.flux 

                idx = np.argsort(bbwave)
                bbwave = np.flip(bbwave[idx])
                bbsflux = np.flip(bbsflux[idx])

                bbflux = integrate.simps(bbwave, bbsflux)

                distcm = self.fit.distance[0]*3.086e+24
                lum = 4*np.pi*(distcm)**2*bbflux

                print('\n\nBLACKBODY MODEL VALIDATION\n'+'#'*80+'\n')
                print('temperature [K]:',str(int(np.round(params[1], decimals=-1))))
                print('peak flux [erg/s/cm2/angstrom]:',np.max(sp.flux))
                print('integrated flux [erg/s/cm2]:',bbflux)
                print('distance [cm]:',distcm)
                print('luminosity [erg/s]:',lum)
                print('luminosity [log(Lsol)]:',np.log10(lum/3.839e33))
                print('\n\n')

            ax.plot(sp.wave, sp.flux, label=title, zorder=1,
                color=palette[i+2], linewidth=0.6*self.figsize)

            if model=='rsg' and not kwargs['simple_rsg']:
                if len(params)==2:
                    params = [0.01,params[0],params[1],200.]
                else:
                    dust_temp = int(truncate(params[3],-2))

                sp1 = self.fit.create_rsg(*params, sptype='star')
                sp2 = self.fit.create_rsg(*params, sptype='dust')
                sp3 = self.fit.create_rsg(*params, sptype='intrinsic')

                if fluxspace: sp1 = self.convert_to_fluxspace(sp1, ext)
                if fluxspace: sp2 = self.convert_to_fluxspace(sp2, ext)
                if fluxspace: sp3 = self.convert_to_fluxspace(sp3, ext)

                ax.plot(sp1.wave, sp1.flux, label='Reddened RSG', zorder=2,
                    color='blue', linewidth=0.6*self.figsize)

                dust_label = '{0} K dust'
                dust_temp = int(truncate(params[3],-1))
                dust_label = dust_label.format(dust_temp)
                ax.plot(sp2.wave, sp2.flux, label=dust_label, zorder=1,
                    color='red', linewidth=0.6*self.figsize)

                ax.plot(sp3.wave, sp3.flux, label='Unreddened RSG', zorder=0,
                    color='plum', linewidth=0.6*self.figsize)

        ylim = [0.5*np.min(flux-fluxerr), 1.2*np.max(flux+fluxerr)]
        xlim = [0.7*np.min(wave-width), 1.2*np.max(wave+width)]
        ax.set_yscale('log')
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        ylim_range = np.arange(int(np.ceil(np.log10(ylim[0]))),
                               int(np.floor(np.log10(ylim[1])))+1)

        if len(ylim_range)==1:
            ylim_range = list(ylim_range)+\
                [ylim_range[0]+0.30102999566]+\
                [ylim_range[0]+0.47712125472]+\
                [ylim_range[0]+0.60205999133]

        maxflux_power = int(np.round(np.log10(1.4*np.max(flux+fluxerr))))
        ytitle = self.create_flux_title(maxflux_power)
        self.setup_axis_titles(ax, 'wavelength', ytitle)

        # Set ticklabels to reflect updated y-axis scale
        ytickval = []
        yticklabel = []
        for val in ylim_range:
            if 10.**val > ylim[1] or 10.**val < ylim[0]: continue
            ytickval.append(10.**val)
            labelval=float('%7.6f'%(10.**val/(10.0**maxflux_power)))
            yticklabel.append(str(labelval))
        ax.set_yticks(ytickval)
        ax.set_yticklabels(yticklabel)

        ytickval = []
        yticklabel = []
        for val in ax.get_yticks(minor=True):
            if val > ylim[1] or val < ylim[0]: continue
            ytickval.append(val)
            labelval=float('%7.6f'%(val/(10.0**maxflux_power)))
            yticklabel.append('')
        ax.set_yticks(ytickval, minor=True)
        ax.set_yticklabels(yticklabel, minor=True)

        if xlim[1] > 3.0e4:
            if xlim[1]>1.0e5: xlim[1]=1.0e5
            ax.set_xlim(xlim)
            ax.set_xscale('log')

            ax.set_xlabel('Observer-frame Wavelength ('+r'$\mu$m'+')')
            ax.set_xticks([3000,10000,30000])
            ax.set_xticklabels(['0.3','1.0','3.0'])

        legend = ax.legend(loc='best', fontsize=3.6*self.figsize)

        # Outfile
        outfile = self.options.objname.strip() + '_sed.eps'
        outfile = os.path.join(self.outdir, outfile)
        self.close_plot(outfile)

    def load_yoon(self, file):
        header = ('model','mass','log_L','reff','log_Teff','helium','co','henv',
            'mH','mHe','ysf','mdot_tr','mdot_wind','frl','SN')
        table = Table.read(file, format='ascii')
        for i,h in enumerate(header):
            table.rename_column('col{0}'.format(i+1), h)
        return(table)

    def load_tracks(self, model_type, params, feh=0.0, clip_early=0.0,
        clip_tail=4.8):
        # Get feh key
        feh_str = str(int(feh*100)).zfill(4)

        all_models = None
        if model_type=='mist':
            feh_str = str(int(feh*100)).zfill(4)

            for m in params:

                mass_str = str(int(m*10000)).zfill(7)

                # We just want theoretical data, so use arbitrary inst/det
                directory = self.fit.dirs['mist']+'FEH_{0}/WFC3/'.format(feh_str)
                fullfile = directory + mass_str + self.fit.files['mist']['suffix']['cmd']

                if not os.path.exists(fullfile):
                    warning = 'WARNING: data do not exist for:{0}'
                    print(warning.format(directory))
                    continue

                cmd_table = ascii.read(fullfile, header_start=14)
                model_table = cmd_table['star_age','log_Teff','log_L']
                model_table.add_column(Column([m]*len(model_table),name='mass'))

                mask = model_table['star_age'] > clip_early
                model_table = model_table[mask]

                mask = model_table['log_Teff'] < clip_tail
                model_table = model_table[mask]

                if not all_models:
                    all_models = model_table
                else:
                    all_models = vstack([all_models,model_table])

        elif model_type=='bpass':
            feh = 0.014
            feh_str = 'z'+str(int(feh*1000)).zfill(3)
            basedir = self.fit.dirs['bpass']+feh_str+'/'

            for model in params:
                modelstr = 'sneplot-{0}-{1}-{2}-{3}'
                if model[2]%1==0: fmtcode='%7.0f'
                else: fmtcode='%7.1f'
                modelstr = modelstr.format(feh_str,int(model[0]),
                    '%7.1f'%model[1],fmtcode%model[2])

                modelfile = basedir + modelstr
                modelfile = modelfile.replace(' ','')

                if not os.path.exists(modelfile):
                    print('{0} model does not exist!'.format(modelfile))
                    return(None)

                all_models = Table.read(modelfile, format='ascii')
                mass = Column([model[0]]*len(all_models),name=('mass'))

                all_models = all_models['col4','col5','col6','col38','col17',
                    'col40','col44']
                all_models.rename_column('col5','log_L')
                all_models.rename_column('col4','log_Teff')
                all_models.rename_column('col6','star1_mass')
                all_models.rename_column('col17','H_mass')
                all_models.rename_column('col38','star2_mass')
                all_models.rename_column('col40','mass_loss_rate')
                all_models.rename_column('col44','rlof_rate')

                all_models.add_column(mass)

        return(all_models)

    def load_progenitors(self, pfile):
        header = ('name','type','log_T','e_log_T', 'log_L','e_log_L')
        data = Table.read(pfile, names=header, format='ascii')
        for i,row in enumerate(data):
            if row['type']=='II-P' or row['type']=='II-L':
                data[i]['type']='II'
        return(data)

    def hr(self, mode='single', progenitors='progenitors.dat',
        add_progenitors=True,add_yoon=False,add_lbvs=False,add_data=True,
        add_rsg=False,limits=False, **kwargs):

        fig, ax = self.setup_plot()
        self.setup_axis_titles(ax, 'log_T', 'log_L')

        # Get range of model track data
        if 'single' in mode:
            tracks = self.load_tracks('mist', [8,10,17,23,30,40,60],
                clip_early=3.0e4, clip_tail=4.67)
        elif 'rsg' in mode:
            tracks = self.load_tracks('mist', [8,10,12,15,18,25],
                clip_early=3.0e4, clip_tail=4.67)
        elif 'binary' in mode:
            tracks = self.load_tracks('bpass', [[19,0.1,0.6]])

        ylim = [0.92*np.min(tracks['log_L'].data),
                1.02*np.max(tracks['log_L'].data)]
        xlim = [1.045*np.max(tracks['log_Teff'].data),
                0.975*np.min(tracks['log_Teff'].data)]

        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        self.setup_ticks(ax)

        if add_rsg:
            rect = mpatches.Rectangle((3.45, 4.0), 0.3, 1.6, linewidth=3,
                edgecolor=(0,0,0,0.5), facecolor=(1,0.8,0.82,0.5))
            rect.zorder = 0
            ax.add_artist(rect)

        for i,mass in enumerate(sorted(np.unique(tracks['mass']))):

            color = palette[i%len(palette)]
            color=black

            mass_track = tracks[tracks['mass']==mass]

            # Get position where track is the hottest to add label
            idx = np.argmax(mass_track['log_Teff'].data)
            label_teff = mass_track['log_Teff'].data[idx]
            label_logL = mass_track['log_L'].data[idx]

            ax.plot(mass_track['log_Teff'].data, mass_track['log_L'].data,
                color=color, linewidth=0.4*self.figsize,zorder=2)

            if mode=='binary':
                terminal = mass_track[-1]

                term_pos = [terminal['log_Teff'],terminal['log_L']]

                ax.plot(term_pos[0],term_pos[1],marker='*',
                    color=color,ms=6*self.figsize)

                ax.text(term_pos[0]-0.0035*self.figsize*(xlim[1]-xlim[0]),
                    term_pos[1]+0.0040*self.figsize*(ylim[1]-ylim[0]), 'SN',
                    horizontalalignment='center', verticalalignment='center',
                    color=color)

                mask = mass_track['H_mass']<1.2*np.min(mass_track['H_mass'])
                h_poor = mass_track[mask]

                h_free_idx = np.argmax(h_poor['H_mass'].data)
                h_free = h_poor[h_free_idx]

                h_free_pos = [h_free['log_Teff'],h_free['log_L']]

                rlof_on_idx = np.argmin(mass_track['rlof_rate'].data)

                idx = np.argmin(mass_track['log_Teff'])
                rlof_on = mass_track[rlof_on_idx]

                rlof_on_pos = [rlof_on['log_Teff'],rlof_on['log_L']]

                ax.plot(rlof_on_pos[0],rlof_on_pos[1],marker='s',
                    color=color,ms=4*self.figsize)

                ax.text(rlof_on_pos[0]+0.0185*self.figsize*(xlim[1]-xlim[0]),
                    rlof_on_pos[1]+0.0070*self.figsize*(ylim[1]-ylim[0]),
                    r'$\mathrm{RLOF}>10^{-5}~M_{\odot}~\mathrm{yr}^{-1}$',
                    horizontalalignment='center', verticalalignment='center',
                    color=color)

            if 'single' in mode:
                idx = np.argmax(mass_track['log_Teff'].data)
                label_teff = mass_track['log_Teff'].data[idx]
                label_logL = mass_track['log_L'].data[idx]
                textpos = [label_teff+0.0065*self.figsize*(xlim[0]-xlim[1]),
                    label_logL-0.004*self.figsize*(ylim[1]-ylim[0])]
            elif 'rsg' in mode:
                idx = -1
                label_teff = mass_track['log_Teff'].data[idx]
                label_logL = mass_track['log_L'].data[idx]
                textpos = [label_teff-0.0030*self.figsize*(xlim[0]-xlim[1]),
                    label_logL+0.0011*self.figsize*(ylim[1]-ylim[0])]
            elif 'binary' in mode:
                idx = 0
                label_teff = mass_track['log_Teff'].data[idx]
                label_logL = mass_track['log_L'].data[idx]
                textpos = [label_teff+0.0065*self.figsize*(xlim[0]-xlim[1]),
                    label_logL-0.004*self.figsize*(ylim[1]-ylim[0])]

                ylim = [0.980*np.min(mass_track['log_L'].data),
                        1.060*np.max(mass_track['log_L'].data)]
                xlim = [1.020*np.max(mass_track['log_Teff'].data),
                        0.980*np.min(mass_track['log_Teff'].data)]

                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

            ax.text(textpos[0], textpos[1], str(mass)+r'$~M_{\odot}$',
                color=color, zorder=5,
                horizontalalignment='center', verticalalignment='center',
                bbox=dict(facecolor='white', edgecolor=black))

        if add_progenitors:
            progdata = self.load_progenitors(self.fit.dirs['data']+progenitors)

            plot_types = [{'type':['II'],'color':red,'marker':'s','name':'SN II'},
                {'type':['IIb'],'color':green,'marker':'o','name':'SN IIb'},
                {'type':['Ib','Ic'],'color':blue,'marker':'D','name':'SN Ib/c'}]

            for ptype in plot_types:
                ax.errorbar([],[],marker=ptype['marker'],
                    ms=2*self.figsize,color=ptype['color'],
                    linewidth=0.4*self.figsize,markeredgecolor=black,zorder=10,
                    markeredgewidth=0.4*self.figsize,label=ptype['name'])


            for ptype in plot_types:
                for row in progdata:
                    if row['type'] in ptype['type']:
                        if row['e_log_L']==0.0: 
                            uplims=[1]
                            lum_err = 0.1
                            # Note - skipping limits
                            continue
                        else: 
                            uplims=[0]
                            lum_err=row['e_log_L']

                        ax.errorbar([row['log_T']],[row['log_L']],
                            xerr=[row['e_log_T']],yerr=[lum_err],
                            uplims=uplims, color=ptype['color'],
                            linewidth=0.4*self.figsize,marker=ptype['marker'],
                            ms=2*self.figsize,zorder=1,capsize=0.5*self.figsize,
                            markeredgecolor=black,
                            markeredgewidth=0.4*self.figsize)

            if 'rsg' in mode:
                legend = ax.legend(loc='lower right',fontsize=4.8*self.figsize)
            else:
                legend = ax.legend(loc='lower left',fontsize=4.8*self.figsize)

        if 'rsg' in mode:
            self.fit.model_type = 'rsg'
            self.load_sed(self.fit.model_type, self.fit.filename)
            reader = self.fit.load_backend(self.fit.model_type, self.fit.phottable)
            sample = np.array(reader.get_chain(flat=True))
            prob = np.array(reader.get_log_prob(flat=True))
            blob = np.array(reader.get_blobs(flat=True))

            params = self.fit.model_fit_params
            ndim = len(params)

            logL=self.fit.calculate_param_best_fit(sample[:,1], prob, ndim, 
                params[1], show=False)
            Teff=self.fit.calculate_param_best_fit(sample[:,2], prob, ndim, 
                params[2], show=False)

            logL_unc = self.fit.calculate_param_best_fit(sample[:,1], prob, 
                ndim, params[1], return_uncertainty=True, show=False)
            Teff_unc = self.fit.calculate_param_best_fit(sample[:,2], prob, 
                ndim, params[1], return_uncertainty=True, show=False)

            logL_unc_dist = 2.17 * self.fit.distance[1]/self.fit.distance[0]

            xerr = np.array([[np.log10(Teff+Teff_unc[0])-np.log10(Teff),
                    np.log10(Teff)-np.log10(Teff-Teff_unc[0])]])
            yerr = np.array([[logL_unc[0] + logL_unc_dist, 
                             logL_unc[1] + logL_unc_dist]])

            ax.errorbar([np.log10(Teff)],[logL],
                xerr=xerr.transpose(),yerr=yerr.transpose(),color=red,
                linewidth=0.4*self.figsize,marker='*',
                markeredgecolor=goldenrod,
                ms=12*self.figsize,zorder=5,capsize=0.5*self.figsize,
                markeredgewidth=0.4*self.figsize)

        if add_yoon:
            yoondata = self.load_yoon('data/yoon_0.02.dat')

            iibflag = False ; ibflag = False
            for row in yoondata:
                label=None
                color=None
                if 'IIb' in row['SN']:
                    if not iibflag: label='Yoon et al. 2017 IIb'
                    color=green
                    iibflag = True
                elif 'Ib' in row['SN']:
                    if not ibflag: label='Yoon et al. 2017 Ib'
                    color=blue
                    ibflag = True
                else:
                    continue
                ax.plot(row['log_Teff'],row['log_L'],'*',
                    ms=4*self.figsize, label=label, color=color)

            if 'rsg' in mode:
                legend = ax.legend(loc='upper left',fontsize=3.2*self.figsize)
            else:
                legend = ax.legend(loc='upper left',fontsize=3.2*self.figsize)
            ylim[0]=4.45

        if add_lbvs:
            ax.plot([np.log10(14000)],[5.55],'o',color=black)
            ax.plot([np.log10(8750)], [5.51],'o',color=black)
            ax.plot([np.log10(7900)], [5.49],'o',color=black)
            ax.plot([np.log10(7500)], [5.45],'o',color=black)
            ax.plot([np.log10(13800)],[5.42],'o',color=black)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if add_data:
            if limits:
                if self.fit.model_type=='blackbody':
                    # Easiest - Then samples represent (luminosity, temperature)
                    # Generate a contour for the region encompassed by inf val
                    photfile = self.fit.dirs['input']+'2020fqv.dat'

                    self.fit.set_model_type('blackbody', extinction=False)
                    self.fit.phottable = self.fit.import_phottable(photfile)

                    mjd, inst_filt, mag, magerr = self.fit.get_fit_parameters(
                        self.fit.phottable)

                    # These values are hard-coded from 2020fqv - should find a
                    # way to do this automatically
                    ext=(1.67,3.22)
                    self.fit.load_extinction(val=ext)

                    self.fit.load_models(self.fit.phottable)
                    self.fit.backend = self.fit.load_backend('blackbody',
                        self.fit.phottable)

                    def ccolor(r,g,b):
                        return((r/255.,g/255.,b/255., 1.0))

                    ncols=60
                    colarray=np.array([ccolor(l,l,l) for l in
                        np.flip(np.arange(255-ncols,255))])
                    cm = ListedColormap(colarray)

                    gsize=(400,500)

                    grid=np.zeros((gsize[0], gsize[1]))
                    lum_space = np.linspace(ylim[0], ylim[1], gsize[1])
                    tem_space = np.linspace(xlim[0], xlim[1], gsize[0])
                    for ii,temp in enumerate(tem_space):
                        for jj,lum in enumerate(lum_space):
                            model,Av,Rv = self.fit.compute_model_mag(inst_filt,
                                (lum, 10**temp), extinction=ext)
                            if any([mmag>nmag for mmag,nmag in zip(mag,model)]):
                                grid[ii,jj]=1

                    vals=np.where(grid==1)
                    idx=np.argmin(vals[1])

                    levels = np.linspace(0, 1, ncols)

                    ax.contourf(tem_space,lum_space,np.transpose(grid),cmap=cm)
                    ax.contour(tem_space,lum_space,np.transpose(grid),
                        linewidths=6, levels=levels, zorder=5, colors=('red'))

            else:
                data=[np.log10(6800),5.3,0.025,0.2]

                if data[0]+data[2]>xlim[0]:
                    xlim[0] = 1.045*(data[0]+data[2])
                if data[0]-data[2]<xlim[1]:
                    xlim[1] = 0.95*(data[0]-data[2])
                if data[1]+data[3]>ylim[1]:
                    ylim[1] = 1.01*(data[1]+data[3])
                if data[1]-data[3]<ylim[0]:
                    ylim[0] = 0.95*(data[1]-data[3])

                if 'rsg' in mode:
                    xlim[0] = 3.89 ; xlim[1] = 3.4
                    ylim[0] = 3.5 ; ylim[1] = 5.6


                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

                xposscale=0.002
                if mode=='single': xposscale=0.002
                elif mode=='binary': xposscale=0.007

                textpos = [data[0]+xposscale*self.figsize*(xlim[0]-xlim[1]),
                           data[1]+0.0082*self.figsize*(ylim[1]-ylim[0])]



        outfile = self.options.objname.strip() + '_hr-{0}.eps'
        outfile = outfile.format(mode[0])
        outfile = os.path.join(self.outdir, outfile)

        self.close_plot(outfile)

    def cmd(self, bands, mode='single', progenitors='progenitors.dat',
        add_progenitors=True,add_yoon=False,add_lbvs=False,add_data=True,
        add_rsg=False,limits=False, **kwargs):

        fig, ax = self.setup_plot()
        self.setup_axis_titles(ax, 'log_T', 'log_L')

        # Get range of model track data
        if mode=='single':
            tracks = self.load_tracks('mist', [8,10,17,23,30,40,60],
                clip_early=3.0e4, clip_tail=4.67)
        elif mode=='binary':
            tracks = self.load_tracks('bpass', [[19,0.1,0.6]])

        abs_data = tracks[bands[0]].data
        col_data = tracks[bands[0]]-tracks[bands[1]].data

        abs_range = np.max(abs_data)-np.min(abs_data)
        col_range = np.max(col_data)-np.min(col_data)

        ylim = [np.min(abs_data)-0.025*abs_range,
                np.max(abs_data)+0.025*abs_range,]
        xlim = [np.min(col_data)-0.025*col_range,
                np.max(col_data)+0.025*col_range]

        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        self.setup_ticks(ax)

        if add_rsg:
            rect = mpatches.Rectangle((3.45, 4.0), 0.3, 1.6, linewidth=3,
                edgecolor=(0,0,0,0.5), facecolor=(1,0.8,0.82,0.5))
            rect.zorder = 0
            ax.add_artist(rect)

        for i,mass in enumerate(sorted(np.unique(tracks['mass']))):

            color = palette[i%len(palette)]
            color=black

            mass_track = tracks[tracks['mass']==mass]

            # Get position where track is the hottest to add label
            idx = np.argmax(mass_track['log_Teff'].data)
            label_teff = mass_track['log_Teff'].data[idx]
            label_logL = mass_track['log_L'].data[idx]

            ax.plot(mass_track['log_Teff'].data, mass_track['log_L'].data,
                color=color, linewidth=0.4*self.figsize,zorder=1)

            if mode=='binary':
                terminal = mass_track[-1]

                term_pos = [terminal['log_Teff'],terminal['log_L']]

                ax.plot(term_pos[0],term_pos[1],marker='*',
                    color=color,ms=6*self.figsize)
                ax.text(term_pos[0]-0.0035*self.figsize*(xlim[1]-xlim[0]),
                    term_pos[1]+0.0040*self.figsize*(ylim[1]-ylim[0]), 'SN',
                    horizontalalignment='center', verticalalignment='center',
                    color=color)

                mask = mass_track['H_mass']<1.2*np.min(mass_track['H_mass'])
                h_poor = mass_track[mask]

                h_free_idx = np.argmax(h_poor['H_mass'].data)
                h_free = h_poor[h_free_idx]

                h_free_pos = [h_free['log_Teff'],h_free['log_L']]

                rlof_on_idx = np.argmin(mass_track['rlof_rate'].data)

                idx = np.argmin(mass_track['log_Teff'])
                rlof_on = mass_track[rlof_on_idx]

                rlof_on_pos = [rlof_on['log_Teff'],rlof_on['log_L']]

                ax.plot(rlof_on_pos[0],rlof_on_pos[1],marker='s',
                    color=color,ms=4*self.figsize)

                ax.text(rlof_on_pos[0]+0.0185*self.figsize*(xlim[1]-xlim[0]),
                    rlof_on_pos[1]+0.0070*self.figsize*(ylim[1]-ylim[0]),
                    r'$\mathrm{RLOF}>10^{-5}~M_{\odot}~\mathrm{yr}^{-1}$',
                    horizontalalignment='center', verticalalignment='center',
                    color=color)

            if mode=='single':
                idx = np.argmax(mass_track['log_Teff'].data)
                label_teff = mass_track['log_Teff'].data[idx]
                label_logL = mass_track['log_L'].data[idx]
                textpos = [label_teff+0.0065*self.figsize*(xlim[0]-xlim[1]),
                    label_logL-0.004*self.figsize*(ylim[1]-ylim[0])]
            elif mode=='binary':
                idx = 0
                label_teff = mass_track['log_Teff'].data[idx]
                label_logL = mass_track['log_L'].data[idx]
                textpos = [label_teff+0.0065*self.figsize*(xlim[0]-xlim[1]),
                    label_logL-0.004*self.figsize*(ylim[1]-ylim[0])]

                ylim = [0.980*np.min(mass_track['log_L'].data),
                        1.060*np.max(mass_track['log_L'].data)]
                xlim = [1.020*np.max(mass_track['log_Teff'].data),
                        0.980*np.min(mass_track['log_Teff'].data)]

                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

            ax.text(textpos[0], textpos[1], str(mass)+r'$~M_{\odot}$',
                color=color, zorder=1,
                horizontalalignment='center', verticalalignment='center')

        if add_progenitors:
            progdata = self.load_progenitors(self.fit.dirs['data']+progenitors)

            plot_types = [{'type':['II'],'color':red,'marker':'s','name':'SN II'},
                {'type':['IIb'],'color':green,'marker':'o','name':'SN IIb'},
                {'type':['Ib','Ic'],'color':blue,'marker':'D','name':'SN Ib/c'}]

            for ptype in plot_types:
                ax.errorbar([],[],marker=ptype['marker'],
                    ms=2*self.figsize,color=ptype['color'],
                    linewidth=0.4*self.figsize,markeredgecolor=black,zorder=10,
                    markeredgewidth=0.4*self.figsize,label=ptype['name'])


            for ptype in plot_types:
                for row in progdata:
                    if row['type'] in ptype['type']:
                        if row['e_log_L']==0.0: uplims=[1]
                        else: uplims=[0]

                        ax.errorbar([row['log_T']],[row['log_L']],
                            xerr=[row['e_log_T']],yerr=[row['e_log_L']],
                            uplims=uplims, color=ptype['color'],
                            linewidth=0.4*self.figsize,marker=ptype['marker'],
                            ms=2*self.figsize,zorder=5,capsize=0.5*self.figsize,
                            markeredgecolor=black,
                            markeredgewidth=0.4*self.figsize)

            legend = ax.legend(loc='lower left',fontsize=3.2*self.figsize)

        if add_yoon:
            yoondata = self.load_yoon('data/yoon_0.02.dat')

            iibflag = False ; ibflag = False
            for row in yoondata:
                label=None
                color=None
                if 'IIb' in row['SN']:
                    if not iibflag: label='Yoon et al. 2017 IIb'
                    color=green
                    iibflag = True
                elif 'Ib' in row['SN']:
                    if not ibflag: label='Yoon et al. 2017 Ib'
                    color=blue
                    ibflag = True
                else:
                    continue
                ax.plot(row['log_Teff'],row['log_L'],'*',
                    ms=4*self.figsize, label=label, color=color)

            legend = ax.legend(loc='upper left',fontsize=3.2*self.figsize)
            ylim[0]=4.45

        if add_lbvs:
            ax.plot([np.log10(14000)],[5.55],'o',color=black)
            ax.plot([np.log10(8750)], [5.51],'o',color=black)
            ax.plot([np.log10(7900)], [5.49],'o',color=black)
            ax.plot([np.log10(7500)], [5.45],'o',color=black)
            ax.plot([np.log10(13800)],[5.42],'o',color=black)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if add_data:
            if limits:
                if self.fit.model_type=='blackbody':
                    # Easiest - Then samples represent (luminosity, temperature)
                    # Generate a contour for the region encompassed by inf val
                    photfile = self.fit.dirs['input']+'2020fqv.dat'

                    self.fit.set_model_type('blackbody', extinction=False)
                    self.fit.phottable = self.fit.import_phottable(photfile)

                    mjd, inst_filt, mag, magerr = self.fit.get_fit_parameters(
                        self.fit.phottable)

                    # These values are hard-coded from 2020fqv - should find a
                    # way to do this automatically
                    ext=(1.67,3.22)
                    self.fit.load_extinction(val=ext)

                    self.fit.load_models(self.fit.phottable)
                    self.fit.backend = self.fit.load_backend('blackbody',
                        self.fit.phottable)

                    def ccolor(r,g,b):
                        return((r/255.,g/255.,b/255., 1.0))

                    ncols=60
                    colarray=np.array([ccolor(l,l,l) for l in
                        np.flip(np.arange(255-ncols,255))])
                    cm = ListedColormap(colarray)

                    gsize=(400,500)

                    grid=np.zeros((gsize[0], gsize[1]))
                    lum_space = np.linspace(ylim[0], ylim[1], gsize[1])
                    tem_space = np.linspace(xlim[0], xlim[1], gsize[0])
                    for ii,temp in enumerate(tem_space):
                        for jj,lum in enumerate(lum_space):
                            model,Av,Rv = self.fit.compute_model_mag(inst_filt,
                                (lum, 10**temp), extinction=ext)
                            if any([mmag>nmag for mmag,nmag in zip(mag,model)]):
                                grid[ii,jj]=1

                    vals=np.where(grid==1)
                    idx=np.argmin(vals[1])

                    levels = np.linspace(0, 1, ncols)

                    ax.contourf(tem_space,lum_space,np.transpose(grid),cmap=cm)
                    ax.contour(tem_space,lum_space,np.transpose(grid),
                        linewidths=6, levels=levels, zorder=5, colors=('red'))

            else:
                data=[np.log10(6800),5.3,0.025,0.2]

                if data[0]+data[2]>xlim[0]:
                    xlim[0] = 1.045*(data[0]+data[2])
                if data[0]-data[2]<xlim[1]:
                    xlim[1] = 0.95*(data[0]-data[2])
                if data[1]+data[3]>ylim[1]:
                    ylim[1] = 1.01*(data[1]+data[3])
                if data[1]-data[3]<ylim[0]:
                    ylim[0] = 0.95*(data[1]-data[3])

                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

                ax.errorbar([data[0]],[data[1]],xerr=[data[2]],yerr=[data[3]],
                    marker='*', ms=4*self.figsize,color=blue,
                    markeredgecolor=goldenrod,zorder=10,
                    capsize=0.5*self.figsize,markeredgewidth=0.3*self.figsize,
                    linewidth=0.4*self.figsize)

                xposscale=0.002
                if mode=='single': xposscale=0.002
                elif mode=='binary': xposscale=0.007

                textpos = [data[0]+xposscale*self.figsize*(xlim[0]-xlim[1]),
                           data[1]+0.0082*self.figsize*(ylim[1]-ylim[0])]


        self.close_plot('cmd-'+mode+'.eps')


    def setup_plot(self, size=[1.8, 1.5]):

        fig, ax = plt.subplots()
        for i in ax.spines.keys(): ax.spines[i].set_linewidth(0.6*self.figsize)
        fig.set_size_inches(1.8*self.figsize, 1.5*self.figsize)

        return(fig, ax)

    def close_plot(self, name):

        if not name.endswith('.eps'): name = name+'.eps'
        name = name.replace('.eps','.png')
        print('Saving figure:',name)

        plt.tight_layout(pad=self.pad)
        plt.savefig(name, format='png')
        plt.clf()
        plt.close('all')

    def setup_ticks(self, ax):

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.tick_params(direction='in', length=2*self.figsize,
            width=0.6*self.figsize, which='major', axis='both', colors=black,
            pad=self.figsize, top=True, bottom=True, left=True, right=True)
        ax.tick_params(direction='in', length=self.figsize,
            width=0.6*self.figsize, which='minor', axis='both', colors=black,
            pad=0.4*self.figsize, top=True, bottom=True, left=True, right=True)

    def setup_axis_titles(self, ax, xname, yname):

        if xname in self.axis_titles.keys(): xname = self.axis_titles[xname]
        if yname in self.axis_titles.keys(): yname = self.axis_titles[yname]

        ax.set_xlabel(xname, labelpad=self.pad)
        ax.set_ylabel(yname, labelpad=self.pad)

    def lightcurve(self, sed, plot_title='', **kwargs):

        phot = sed.phottable
        magunit = sed.magsystem.replace('mag','').upper().strip()
        ytitle = 'Apparent Magnitude ({0} mag)'.format(magunit)

        if 'discovery_date' in phot.meta.keys():
            xtitle = 'time'
            relmjd = Time(phot.meta['discovery_date']).mjd
        else:
            xtitle = 'time_mjd'
            relmjd = 0.0

        fig, ax = self.setup_plot()
        self.setup_axis_titles(ax, xtitle, ytitle)
        self.setup_ticks(ax)

        ylim = [1.005*np.max(mag+err),  0.995*np.min(mag-err)]
        xlim = [1.01*np.min(date), 0.99*np.max(date)]

        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        for i,filt in enumerate(sorted(np.unique(phot['filter']))):

            # Get data set
            mask = phot['filter']==filt
            date, inst_filt, mag, err = sed.get_fit_parameters(phot[mask])
            color = color_palette[inst_filt[0]]

            # Plot data
            ax.errorbar(date-relmjd, mag, yerr=err, capsize=self.figsize,
                linewidth=0.8*self.figsize, label=filt.upper(), color=color,
                capthick=0.8*self.figsize, zorder=10)

        legend = ax.legend(loc='lower right', fontsize=3.2*self.figsize)

        if not plot_title: plot_title = 'lightcurve'

        self.close_plot(plot_title)

    def plot_corner(self, models, plot_blobs=False, notau=False, **kwargs):
        phot = self.fit.phottable
        dum, inst_filt, mag, magerr = self.fit.get_fit_parameters(phot)

        # Get flux for observations in Janskies
        flux = 3631 * 10**(-0.4 * np.array(mag)) * 1.0e-23
        fluxerr = flux * np.array(magerr)/1.086
        fluxcorr = False

        bps = [self.fit.inst_filt_to_bandpass(i) for i in inst_filt]
        bandpasses = np.array(bps)

        for i,model in enumerate(models):

            print('Plotting model:',model)
            self.load_sed(model, self.fit.filename, notau=notau)

            npad = 0.006

            rc('font',**{'family':'serif','serif':['Times'],
                'size':3.1*self.figsize})
            rc('text',**{'usetex':True, 'color':black})

            self.fit.set_model_type(model, extinction=True)
            self.fit.backend = self.fit.load_backend(model, self.fit.phottable)

            ndim = len(self.fit.model_fit_params)
            samples = np.array(self.fit.backend.get_chain(flat=True))
            prob = np.array(self.fit.backend.get_log_prob(flat=True))

            samples, prob = self.fit.sample_params(samples, prob, ndim,
                downsample=0.8)

            nsamples = len(samples)
            print(f'Plottinng {nsamples} samples for model {model}')

            if plot_blobs:
                blobs = np.array(self.fit.backend.get_blobs(flat=True))
                blobs = blobs[start_idx:]

                Av_step = 0.1
                Rv_step = 0.08333333333

                adjust = rand(blobs.shape[0]*blobs.shape[1]) - 0.5
                adjust = adjust.reshape(blobs.shape)
                adjust[:,1] = adjust[:,1] * 2*Av_step
                adjust[:,0] = adjust[:,0] * 2*Rv_step

                blobs = blobs + adjust

                data = np.hstack([samples, blobs])
            else:
                data = samples

            # Get rid of nan values in samples and blobs
            mask = ~np.isnan(data)
            mask = np.all(mask, axis=1)
            data = data[mask, :]

            if model=='pickles':

                fig = corner.corner(data, bins=30,
                    labels=(self.axis_titles['log_L'],
                            self.axis_titles['log_T'],
                            r'$A_{V}$ (mag)',r'$R_{V}$'),
                    ms=12, title_kwargs={'fontsize': 3.1*self.figsize},
                    hist_kwargs={'linewidth': 0.4*self.figsize},
                    labelpad=npad*self.figsize)

            elif model=='bpass':

                fig = corner.corner(data, bins=30,
                    labels=(self.axis_titles['initial_mass'],
                            self.axis_titles['period'],
                            self.axis_titles['mass_ratio'],
                            r'$A_{V}$ (mag)',r'$R_{V}$'),
                    ms=12, title_kwargs={'fontsize': 3.1*self.figsize},
                    range=[(15,25),(0.1,0.4),(0.1,0.9),(1.6,4.0),(1.0,6.0)],
                    hist_kwargs={'linewidth': 0.4*self.figsize},
                    labelpad=npad*self.figsize)

            elif model=='rsg':
                # Transform tau_V -> A_V
                data[:,0] = 0.79 * data[:,0]

                data_range = []
                for i in np.arange(4):
                    data_range.append((0.95*np.min(data[:,i]),
                        1.05*np.max(data[:,i])))

                fig = corner.corner(data, bins=20,
                    labels=(self.axis_titles['A_V'],
                            self.axis_titles['log_L'],
                            self.axis_titles['Teff'],
                            self.axis_titles['Tdust']),
                    ms=12, title_kwargs={'fontsize': 3.1*self.figsize},
                    range=data_range,
                    hist_kwargs={'linewidth': 0.4*self.figsize},
                    labelpad=npad*self.figsize)
            else:
                print(f'ERROR: {model} not recognized for corner plot!')
                return(None)

            fig.set_size_inches(1.5*self.figsize, 1.5*self.figsize)
            plt.tight_layout(pad=0.1*self.pad, w_pad=self.pad, h_pad=self.pad)

            outfile = self.options.objname + f'_corner_{model}.pdf'
            outfile = self.outdir + outfile

            fig.savefig(outfile, format='pdf')

            subprocess.call(['pdf2ps', outfile, outfile.replace('.pdf','.eps')])


def main():
    sed = sed_plot()

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

    if not os.path.exists(sed.outdir):
        os.makedirs(sed.outdir)

    photfile = sed.fit.parse_photfile(opt.objname)
    if not photfile or not os.path.exists(photfile):
        inpfile = sys.argv[1]
        print(f'ERROR: input photfile {inpfile} does not exist!  Exiting...')
        sys.exit(1)

    sed.load_sed('blackbody', photfile, notau=sed.options.notau)

    for typ in opt.type.split(','):
        if typ not in sed.plot_types.keys():
            typs = ','.join(list(sed.plot_types.keys()))
            print(f'WARNING: plot type {typ} not supported.')
            print('Plot type must be one of {typs}.  Continuing...')

        models = opt.model.split(',')

        func = sed.plot_types[typ]
        kwargs=dict(opt._get_kwargs())
        func(models, **kwargs)

if __name__ == "__main__":
    main()
