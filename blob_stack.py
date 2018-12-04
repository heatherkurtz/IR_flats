#file:blob_stack.py
#Description: Stacks the MDI files based on blob apearances one full image 
#per appearance. This can only be used if the 512 flag in not used when creating the mdi files.
#Date: Oct. 25, 2018
#author: Heather Kurtz

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


def get_expstar(file,mjd_list):
	hdr = fits.getheader(file,0)
	mjd = hdr['EXPSTART']
	mjd_list.append(mjd) 
	return(mjd_list)


#def read_in_bpix(bpix_file):
#	return()


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
#				bpix_data0.append(file)
#
#
		bpix_date_list=[55285.58, 54285.58, 56285.58 ]
		mjd_list=[]
		date_dict={}
		for file in list_of_files:
			mjd_list=get_expstar(file,mjd_list)

			
		for bpix_date in bpix_date_list:

			date_mask=[bpix_date < mjd_list]
			date_dict[bpix_date]=list_of_files[date_mask]
		print(date_dict)
main()































