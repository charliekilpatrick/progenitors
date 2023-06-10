from analysis import *

sed = sed_fitter()
photfile = sed.dirs['input']+'2020fqv.dat'

sed.set_model_type('blackbody', extinction=False)
sed.phottable = sed.import_phottable(photfile)
mjd, inst_filt, mag, magerr = sed.get_fit_parameters(sed.phottable)

ext=(1.67,3.22)
sed.load_extinction(val=ext)

sed.load_models(sed.phottable)
sed.backend = sed.load_backend('blackbody', sed.phottable)

model, Av, Rv = sed.compute_model_mag(inst_filt, (5.189890611688778, 3300.6884332127615),
    extinction=ext)

print(inst_filt)
print(model)
print(mag)
print(model > mag)
print(any([m>n for m,n in zip(mag, model)]))
