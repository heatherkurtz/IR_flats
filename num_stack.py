#name: num_stack.py
#Description: This creates a mask of all the cosmic rays in the DQ and applies it to the image 
#Date: August 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

filters = ['F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W']

current = os.getcwd()
#base_path = '/user/hkurtz/IR_flats/test_f098'
#shange directory to the base path
base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/nodq'
os.chdir(base_path)

for fil in filters:
	fil_path= os.path.join(base_path,fil)
	os.chdir(fil_path)

	#list of files
	list_of_files= glob.glob('*mdi.fits')
	
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
	n_num = fil + '_num.fits'
	fits.writeto(n_num, im, overwrite=True)
	
