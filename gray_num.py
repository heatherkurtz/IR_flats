#name: gray_num.py
#Description: combines images into stack count the images put into each pixel
#Date: August 12, 2018
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

current = os.getcwd()
#shange directory to the base path
base_path = '/grp/hst/wfc3v/hkurtz/sky_flats'
os.chdir(base_path)

#list of files
list_of_files= glob.glob('*/*mdi.fits')

#gets the file size
hdr = fits.getheader(list_of_files[0], 1)
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']
nf = len(list_of_files)
set_data=fits.getdata(list_of_files[0], 1)
#makes empty array to be filled with data
data_array = np.empty((nf, ny, nx), dtype=float)

#read in the data and the DQs from the .fits for both extensions

for i , f in enumerate(list_of_files):

	data=fits.getdata(f, 1)
	data[data > 0] = 1
	data_array[i, :, :] = data
im = np.nansum(data_array, axis=0)
n_num = 'gray_num.fits'
fits.writeto(n_num, im, overwrite=True)
	
