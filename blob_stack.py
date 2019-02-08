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

filters = ['F160W']#[ 'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W' ]

#'F132N', 'F139M', 'F126N', 'F153M', 'F127M','F128N', 'F164N'
#'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W',

def Stack(data_array_1, signal_arr):
	"""
	combines the 3D array of images into one image for mean, median, sum, signal, and std

	:param data_array_1: Numpy array of images.
    :param signal_arr: Numpy array of signal images.
    """
	image_median = np.nanmedian(data_array_1, axis=0)
	image_mean = np.nanmean(data_array_1, axis=0)
	image_sum = np.nansum(data_array_1, axis=0)
	image_signal = np.nansum(signal_arr, axis=0)
	image_std = np.nanstd(data_array_1, axis=0)
	#image_min = np.nanmin(data_array_1, axis=0)

	return (image_mean, image_median, image_sum, image_std, image_signal)


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
		print('working on ', fil)
		fil_path= os.path.join(base_path,fil)
		os.chdir(fil_path)
		#list of files
		list_of_files= glob.glob('*mdi.fits')
		hdr = fits.getheader(list_of_files[0], 1)
		nx = hdr['NAXIS1']
		ny = hdr['NAXIS2']
		#print(len(list_of_files))

#		for file in list_of_files:
#			hdr = fits.getheader(file,0)
#			mjd = hdr('EXPSTART')
#			if bpix_data[0] < mjd:
#				bpix_data0.append(file)
#			if bpix_data[1] < mjd:
#				bpix_data1.append(file)
#
#
		bpix_date_list = read_in_bpix()
		mjd_list = []
		date_dict = {}
		for file in list_of_files:
			mjd_list = get_expstar(file,mjd_list)

		mjd_arr = np.array(mjd_list)

		files_arr = np.array(list_of_files)
		for bpix_date in bpix_date_list:
			#this masks images with an observation date less than the date of the blob aprearances and stores them in a dictionary
			date_mask = [bpix_date < mjd_arr]
			date_dict[bpix_date] = files_arr[date_mask]
		for key in date_dict:
			# creats 3D array of size of the number of images in dictionary
			nf = len(date_dict[key])
			data_array = np.empty((nf, ny, nx), dtype=float)
			sig_3Darr = np.empty((nf, ny, nx), dtype=float)
			for i, f in enumerate(date_dict[key]):
				#reads in the data
				data=fits.getdata(f, 1)
				hdr = fits.getheader(f, 0)
				norm = hdr['NORM']
				exptime = hdr['EXPTIME']
				data_array[i, :, :] = data
				sig_3Darr[i, :, :] = data * norm * exptime
			S_mean,S_median,S_sum,S_std, S_signal=Stack(data_array, sig_3Darr)
			#writing out the files
			n_mean = fil + '_' + str(key) + '_' + 'blob_mean.fits'
			n_median = fil + '_' + str(key) + '_' + 'blob_median.fits'
			n_sum = fil + '_' + str(key) + '_' + 'blob_sum.fits'
			#n_min = fil + '_' + str(key) + '_' + 'blob_min.fits'
			n_std = fil + '_' + str(key) + '_' + 'blob_std.fits'
			n_signal = fil + '_' + str(key) + '_' + 'blob_signal.fits'
			fits.writeto(n_mean, S_mean,overwrite=True)
			fits.writeto(n_median, S_median,overwrite=True)
			fits.writeto(n_sum, S_sum,overwrite=True)
			fits.writeto(n_std, S_std,overwrite=True)
			#fits.writeto(n_min, S_min,overwrite=True)
			fits.writeto(n_signal, S_signal,overwrite=True)
		
main()































