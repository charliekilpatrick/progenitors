from astropy.table import Table
import pysynphot as S
import numpy as np
import os
import matplotlib.pyplot as plt

data = Table.read('data/pickles.dat',format='ascii')

stars = ['K4I','B5III']
scales = [1.0, 0.3]

sps=[]
for star in stars:

    mask = np.array([star in row['SPTYPE'] for row in data])
    filename = data[mask]['FILENAME'][0]

    fullfile = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles',
            'dat_uvk', filename+'.fits')

    if os.path.exists(fullfile):
            sp = S.FileSpectrum(fullfile)
            sps.append(sp)


newsp = scales[0]*sps[0] + scales[1]*sps[1]
newsp = newsp * 1.0e8 * 1.1
print(newsp.fluxunits)
newsp1 = scales[0]*sps[0] * 1.0e8 * 1.1
newsp2 = scales[1]*sps[1] * 1.0e8 * 1.1

flux = newsp.flux
wave = newsp.wave

plt.plot(wave, flux, color='green')
plt.plot(newsp1.wave, newsp1.flux, color='red')
plt.plot(newsp2.wave, newsp2.flux, color='blue')
plt.ylabel(r'Flux')
plt.xlabel(r'Wavelength ($\AA$)')
plt.xlim(3000,10000)
plt.ylim(0, 1)
plt.savefig('test.png')
