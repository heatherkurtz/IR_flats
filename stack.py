#Filename: mask.py
#Description: This creates a mask of all the cosmic rays in the DQ and applies it to the image 
#Date: August 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

filters = [ 'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W' ]

#'F132N', 'F139M', 'F126N', 'F153M', 'F127M','F128N', 'F164N'
#'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W',

def Stack(data_array_1):
	image_median = np.nanmedian(data_array_1, axis=0)
	image_mean = np.nanmean(data_array_1, axis=0)
	image_sum = np.nansum(data_array_1, axis=0)
	image_std = np.nanstd(data_array_1, axis=0)
	image_min = np.nanmin(data_array_1, axis=0)

	return (image_mean, image_median, image_sum, image_std, image_min)


#def number_stack(data_array_1):
	#get number of good values per pixel




def main():

	#gets current directory
	current = os.getcwd()
	#base_path = '/user/hkurtz/IR_flats/test_f098'
	base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/nodq/'
	#shange directory to the base path
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
		
			data_1=fits.getdata(f, 1)
			data_array[i, :, :] = data_1
		S_mean,S_median,S_sum,S_std,S_min=Stack(data_array)
		n_mean = fil + '_mean.fits'
		n_median = fil + '_median.fits'
		n_sum = fil + '_sum.fits'
		n_min = fil + '_min.fits'
		n_std = fil + '_std.fits'
		fits.writeto(n_mean, S_mean,overwrite=True)
		fits.writeto(n_median, S_median,overwrite=True)
		fits.writeto(n_sum, S_sum,overwrite=True)
		fits.writeto(n_std, S_std,overwrite=True)
		fits.writeto(n_min, S_min,overwrite=True)
	
main()




