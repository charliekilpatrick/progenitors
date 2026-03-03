"""Example: mass-distribution modeling (progenitors.davies18.distribution).
Run from repo root: python examples/example_simulation.py"""
import numpy as np
from progenitors.davies18 import distribution

generate_distribution = distribution.generate_distribution
generate_sample = distribution.generate_sample
calculate_chi2 = distribution.calculate_chi2

# generate_distribution
masses, probs = generate_distribution(10.0, 15.0)
assert len(masses) == int((15 - 10) / 0.1) and len(probs) == len(masses)
assert np.all(probs > 0) and masses[0] == 10.0 and masses[-1] == 15.0
# generate_sample
m = np.array([1.0, 2.0, 3.0])
p = np.array([0.2, 0.5, 0.3])
samples = generate_sample(m, p, 10)
assert len(samples) == 10 and np.all(np.isin(np.array(samples).flatten(), m))
# calculate_chi2
inp = np.array([10.0, 12.0, 14.0])
sim = np.array([10.1, 11.9, 14.2])
lims = np.array([1.0, 1.0, 1.0])
chi2 = calculate_chi2(inp, sim, lims)
assert chi2 >= 0 and np.isfinite(chi2)
assert calculate_chi2(inp, inp, lims) == 0.0
print("progenitors.davies18.distribution OK")
