import numpy as np
import random

def generate_distribution(m_min, m_max):

    size = int((m_max - m_min)/0.1)
    masses = np.linspace(m_min, m_max, size)

    probs = masses**-1.35

    return(masses, probs)

def generate_sample(masses, probs, size):

    output = []
    for i in np.arange(size):
        output.append(random.choices(masses, probs))

    return(output)

def calculate_chi2(input_masses, sim_masses, lims):

    chi2 = np.sum(lims*(input_masses-sim_masses)**2) * 1.0/len(input_masses)**2
    return(chi2)

N=100000
masses, lims = np.loadtxt('input_masses.txt', dtype=float, unpack=True)
print(len(masses))
idx = np.argsort(masses)
masses = masses[idx]
lims = lims[idx]
vals = []
chi = []
for i in np.arange(N):

    print(i)

    m_min = np.random.uniform(6.0, 10.0)
    m_max = np.random.uniform(15.0, 35.0)
    mask = masses < m_max
    inp_masses = masses[mask]
    inp_lims = lims[mask]
    vals.append((m_min, m_max))

    sim_masses, prob = generate_distribution(m_min, m_max)

    samp_masses = generate_sample(sim_masses, prob, len(inp_masses))
    samp_masses = sorted(samp_masses)

    chi.append(calculate_chi2(inp_masses, samp_masses, inp_lims))

chi = np.array(chi)
vals = np.array(vals)
chi = chi/np.min(chi)

for s in [1.0, 2.3, 3.5]:
    mask = chi < 1.0 + s
    mask_vals = vals[mask]
    m_min = np.array([v[0] for v in mask_vals])
    m_max = np.array([v[1] for v in mask_vals])

    print(s)
    print(np.min(m_max))
    print(np.max(m_max))
    print(np.median(m_max))
