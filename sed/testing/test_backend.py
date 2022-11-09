import emcee
import numpy as np

backfile='../data/backends/2020fqv_blackbody.h5'
newname='15495456422528650199174'

backend = emcee.backends.HDFBackend(backfile, name=newname)

flat_samples = np.array(backend.get_chain(flat=True))
flat_prob = np.array(backend.get_log_prob(flat=True))
flat_blobs = np.array(backend.get_blobs(flat=True))

mask=~np.isinf(flat_prob)

print(len(flat_prob[mask]))
