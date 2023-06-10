from astropy.table import Table
import glob
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

for file in glob.glob('data/mist/FEH_0000/ACS_WFC/*.cmd'):

    table = Table.read(file, format='ascii', header_start=14)

    plt.plot(table['ACS_WFC_F555W']-table['ACS_WFC_F814W'],
        table['ACS_WFC_F555W'])

plt.plot([0.3],[-9.1], 'o')

plt.savefig('cmd.png')


