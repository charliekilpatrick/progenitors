import util

import warnings
warnings.filterwarnings('ignore')

row = {'RA': '01:36:48.161', 'Dec': '+15:45:31.00'}

output = util.get_jwst_data(row)
print(output)
