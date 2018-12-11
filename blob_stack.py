#file:blob_stack.py
#Description: Stacks the MDI files based on blob apearances one full image 
#per appearance. This can only be used if the 512 flag in not used when creating the mdi files.
#Date: Oct. 25, 2018
#author: Heather Kurtz

from astropy.io import ascii
from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

filters = [ 'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W' ]

#'F132N', 'F139M', 'F126N', 'F153M', 'F127M','F128N', 'F164N'
#'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W',

def _Stack(data_array_1, signal_arr):
	image_median = np.nanmedian(data_array_1, axis=0)
	image_mean = np.nanmean(data_array_1, axis=0)
	image_sum = np.nansum(data_array_1, axis=0)
	image_signal = np.nansum(signal_arr, axis=0)
	image_std = np.nanstd(data_array_1, axis=0)
	image_min = np.nanmin(data_array_1, axis=0)

	return (image_mean, image_median, image_sum, image_std, image_min, image_signal)


def get_expstar(file,mjd_list):
	hdr = fits.getheader(file,0)
	mjd = hdr['EXPSTART']
	mjd_list.append(mjd) 
	return(mjd_list)


def read_in_bpix():
	filt_table = ascii.read('/user/hkurtz/IR_flats/bpixtab_summary.txt', data_start=2, delimiter=' ')
	usaftermjd = filt_table['useafter_mjd']
	return (usaftermjd)


def main():

	#gets current directory
	current = os.getcwd()
	#base_path = '/user/hkurtz/IR_flats/test_f098'
	base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/no512/'
	#shange directory to the base path
	os.chdir(base_path)
	
	for fil in filters:
		fil_path= os.path.join(base_path,fil)
		os.chdir(fil_path)
		#list of files
		list_of_files= glob.glob('*mdi.fits')

#		for file in list_of_files:
#			hdr = fits.getheader(file,0)
#			mjd = hdr('EXPSTART')
#			if bpix_data[0] < mjd:
#				bpix_data0.append(file)
#			if bpix_data[1] < mjd:
#				bpix_data1.append(file)
#
#
		bpix_date_list=read_in_bpix()
		mjd_list=[]
		date_dict={}
		for file in list_of_files:
			mjd_list=get_expstar(file,mjd_list)

		mjd_arr=np.array(mjd_list)
		files_arr=np.array(list_of_files)
		for bpix_date in bpix_date_list:

			date_mask=[bpix_date < mjd_arr]
			date_dict[bpix_date]=files_arr[date_mask]
		print(date_dict)
main()































