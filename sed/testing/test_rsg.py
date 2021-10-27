import analysis
import pysynphot as S

sed = analysis.sed_fitter(interpolate=False)
photfile = sed.dirs['input']+'2020jfo.dat'

sed.set_model_type('rsg', extinction=False)

sed.phottable = sed.import_phottable(photfile)
mjd, inst_filt, mag, magerr = sed.get_fit_parameters(sed.phottable)
bandpasses = [sed.inst_filt_to_bandpass(i) for i in inst_filt]

sed.load_extinction(val=(0.186, 3.1))
sed.load_models(sed.phottable)

ext = sed.host_ext
print(ext)

theta = (3.276690394596393, 4.680673923430319, 4877.866202492276, 1925.0427941032867)

mags = sed.compute_rsg_mag(inst_filt, *theta)
sp = sed.create_rsg(*theta)
print(mags)
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

bandpasses = sed.get_bandpasses(inst_filt)
mags = sed.compute_pysynphot_mag(sp, bandpasses)
print(mags)
print(mags-mag)

model, Av, Rv = sed.compute_model_mag(inst_filt, theta, extinction=ext)
chi2, Av, Rv = sed.log_likelihood(theta, inst_filt, mag, magerr,
            extinction=ext)

print(inst_filt)
sed.print_model_results(inst_filt, model, mag, magerr, theta, ext)
