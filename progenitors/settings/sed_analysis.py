"""
SED model fit parameter names by model type.

Used by: progenitors.sed.analysis (sed_fitter)
Location: progenitors/settings/sed_analysis.py
"""
# Model type -> list of parameter names to fit
MODEL_FIT_PARAMS = {
    "bpass": ["mass", "ratio", "period"],
    "mist": ["mass", "age"],
    "mist_terminal": ["mass"],
    "blackbody": ["luminosity", "temperature"],
    "pickles": ["luminosity", "temperature"],
    "rsg": ["tau_V", "luminosity", "temperature", "dust_temp"],
    "rsg_sil_2": ["tau_V", "luminosity", "temperature", "dust_temp"],
    "rsg_grp_2": ["tau_V", "luminosity", "temperature", "dust_temp"],
    "rsg_sil_10": ["tau_V", "luminosity", "temperature", "dust_temp"],
    "rsg_grp_10": ["tau_V", "luminosity", "temperature", "dust_temp"],
    "grams": ["luminosity", "temperature", "mlr"],
}
