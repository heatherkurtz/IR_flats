#file: blob_plots.py
#start with plot of blob 25

import matplotlib.pyplot as plt
import csv
import pandas as pd 
import glob
from astropy.io import ascii 
import numpy as np 
import os
from astropy.io import fits
from scipy import ndimage





def main():
	list_filters = ['F160W','F098M','F105W','F125W','F140W','F110W']
	current = os.getcwd()


	for filt in list_filters:
		pflat= '/grp/hst/wfc3c/mack/ir_flat/PIRZKAL_2017/test_delivery/' + filt + '_new_pfl.fits'
		pflat_data = fits.getdata(pflat)
		path = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_no512_run/' + filt +'/'#nov_4_no512_run
		os.chdir(path)
		deltas = glob.glob('*delta_ground_blobs_no_nan_median_delta.fits')
		for f in deltas:
			data = fits.getdata(f)
			data[data == 1] = 0
			#print(data[715:720,570:575])
			#mult_delta = data * pflat_data
			#mult_delta[mult_delta == 0] = 1
			div_delta = data / pflat_data
			div_delta[div_delta == 0] = 1

			#result = ndimage.generic_filter(div_delta, np.nanmedian, size=3, mode='constant', cval=np.NaN)
			#div_delta[np.isnan(div_delta)] = result[np.isnan(div_delta)]

			
			#print(div_delta[715:720,570:575])
			#div_delta[div_delta is np.nan] = 1
			div_name = f[:-5] + "_delta_divid_newpflat_v3.fits"
			#mult_name = f[:-5] + "_mult_ones_test_median.fits"
			fits.writeto(div_name,div_delta,overwrite=True)
			#fits.writeto(mult_name,mult_delta,overwrite=True)
			#mask_nans = 
			#dq = np.copy(div_delta)
			#blob_im1full = fits.getdata('/user/holszewski/IR_flats/ground_blobs.fits')
			#blob_im1 = blob_im1full[5:1019,5:1019]
			#dq[dq !=1] = 512
			#f_dq = dq + blob_im1full
			#f_dq[f_dq < 512] = 0
			#f_dq[f_dq > 512] = 512
			#dq_name = f[:-5] + "_dq_array_test_median.fits"
			#fits.writeto(div_name,div_delta,overwrite=True)
			#fits.writeto(mult_name,mult_delta,overwrite=True)
			#fits.writeto(dq_name,f_dq,overwrite=True)


main()








