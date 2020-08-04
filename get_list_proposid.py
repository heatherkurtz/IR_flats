#file: get_list_proposid.py
import glob
from astropy.io import fits

#path = '/grp/hst/wfc3v/hkurtz/sky_flats/input_data'

files = glob.glob('/grp/hst/wfc3v/hkurtz/sky_flats/MAST_data/*raw.fits')
list_ids = []
for f in files:
	hdr = fits.getheader(f)
	ID = hdr['PROPOSID']
	list_ids.append(ID)


ID_list = list(set(list_ids))
print(ID_list)

