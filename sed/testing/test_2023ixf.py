import analysis
import pysynphot as S
import numpy as np

sed = analysis.sed_fitter(interpolate=False)
photfile = sed.dirs['input']+'2023ixf.dat'

sed.set_model_type('rsg', extinction=False)

sed.phottable = sed.import_phottable(photfile)
mjd, inst_filt, mag, magerr = sed.get_fit_parameters(sed.phottable)
bandpasses = [sed.inst_filt_to_bandpass(i) for i in inst_filt]

sed.load_extinction(val=(0.186, 3.1))
sed.load_models(sed.phottable)

ext = sed.host_ext

for i in np.arange(31):

    tau=i*0.2

    theta = (tau, 4.74, 3920, 880)

    sp = sed.create_rsg(*theta)
    distance = sed.distance[0]
    # Assume distance in Mpc, so D/10pc = distance * 1e5
    distance = distance * 1e5
    # This is for host extinction
    extinction1 = sed.extinction_law(sp.wave, ext[0], ext[1])
    # This is for MW extinction
    Av = sed.mw_ebv * 3.1
    extinction2 = sed.extinction_law(sp.wave, Av, 3.1)
    # Now calculate observed flux and re-generate spectrum
    obsflux = sp.flux*extinction1.flux*extinction2.flux
    obsflux = obsflux/(distance**2)
    sp = S.ArraySpectrum(sp.wave, obsflux, fluxunits='flam')

    with open('2023ixf_tau{0}.dat'.format('%.3f'%tau), 'w') as f:
        f.write('wave flux \n')
        for i in np.arange(len(sp.wave)):
            f.write('{0} {1} \n'.format(
                '%.5f'%sp.wave[i], '%.5e'%sp.flux[i]))
