import numpy as np
import sys
import shutil
import os

filename=sys.argv[1]
if not os.path.exists(filename):
    print(f'ERROR: {filename} does not exist.  Exiting...')
dirname = os.path.dirname(filename)
basefile = os.path.basename(filename)
tmpfile = os.path.join(dirname, 'tmp')

with open(tmpfile, 'w') as f:
    wave, trans = np.loadtxt(filename, unpack=True, dtype=float)

    idx = np.argsort(wave)
    wave = wave[idx]
    trans = trans[idx]

    mask = trans < 0.0
    trans[mask] = 0.0

    if trans[0]!=0.0:
        w0 = wave[0] - (wave[1]-wave[0])
        wave = np.array([w0] + list(wave))
        trans = np.array([0.0] + list(trans))
    if trans[-1]!=0.0:
        w1 = wave[-1] + (wave[-1]-wave[-2])
        wave = np.array(list(wave) + [w1])
        trans = np.array(list(trans) + [0.0])

    vals, idx = np.unique(wave, return_index=True)
    wave = wave[idx]
    trans = trans[idx]

    for w,t in zip(wave,trans):
        line = '%.3f'%w
        line += ' '
        line += '%.5f'%t
        f.write(line+'\n')

shutil.move(tmpfile, filename)
