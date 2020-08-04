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
from bisect import bisect



def read_bpix_file():
    filt_table = ascii.read('/user/holszewski/IR_flats/bpixtab_summary.txt', data_start=1, delimiter=',')
    bpixtab = filt_table['col1']
    usaftermjd = filt_table['col3']
    return (bpixtab,usaftermjd)

def weak_blob_mask(data):
    blob_mask = fits.getdata('/user/holszewski/IR_flats/test_ground_blobs.fits')
    blob_mask[blob_mask > 0] = 1
    data[blob_mask > 0] = np.mean(sigma_clip(data, sigma=5))

def match_date(date,useatermjd,bpixtab): 
    dat = float(date)
    index = bisect(useatermjd, date) -1
    #corin=index -1
    match = bpixtab[index]
    match_name = match
    return(match_name)

def dq_mask(dq):
    bit_mask = (512)
    dq0 = np.bitwise_and(dq, np.zeros(np.shape(dq), 'Int16') + bit_mask)
    dq0[dq0 > 0] = 1
    return dq0




def main():
	list_filters = ['F160W','F098M','F105W','F125W','F140W','F110W']
	current = os.getcwd()


	for filt in list_filters:
		pflat= '/grp/hst/wfc3c/mack/ir_flat/PIRZKAL_2017/test_delivery/' + filt + '_pfl.fits'
		pflat_data = fits.getdata(pflat)
		path = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_no512_run/' + filt +'/'#nov_4_no512_run
		os.chdir(path)
		deltas = glob.glob('*_delta.fits')
		for f in deltas:
			data = fits.getdata(f)
			data[data == 1] = 0
			mult_delta = data * pflat_data
			mult_delta[mult_delta == 0] = 1
			div_delta = data / pflat_data
			div_delta[div_delta == 0] = 1
			for pix in div_delta:
				if pix is np.nan:
					pix=1
			#div_delta[div_delta is np.nan] = 1
			div_name = f[:-5] + "_div_ones.fits"
			mult_name = f[:-5] + "_mult_ones.fits"
			fits.writeto(div_name,div_delta,overwrite=True)
			fits.writeto(mult_name,mult_delta,overwrite=True)
			#mask_nans = 
			dq = np.copy(div_delta)
			bpix, bpix_date_list = read_bpix_file()
			#bpix_dict = {}
			date = f[6:-11]
			print(date)

			bpixtab = match_date(date,bpix_date_list,bpix)
			blob_im_name = '/user/holszewski/IR_flats/' + bpixtab[:-11] + ".fits"
			blob_im_full = fits.getdata(blob_im_name)
			blob_im = blob_im_full[5:1019,5:1019]
			blob_mask = dq_mask(blob_im)
			#blob_im[blob_im < 511] = 0
			#blob_im[blob_im > 511] = 1
			blob_im1full = fits.getdata('/user/holszewski/IR_flats/ground_blobs.fits')
			blob_im1 = blob_im1full[5:1019,5:1019]
			blob_mask1 = dq_mask(blob_im1)
			#blob_im1[blob_im1 < 511] = 0
			#blob_im1[blob_im1 > 511] = 1
			f_blob=blob_mask1+blob_mask
			f_blob[f_blob < 1] = 0
			f_blob[f_blob > 0] = 1


			#dq[dq != 1] = 512
			#dq[dq == 1] = 0
			dq_name = f[:-5] + "_dq_array2.fits"
			#fits.writeto(div_name,div_delta,overwrite=True)
			#fits.writeto(mult_name,mult_delta,overwrite=True)
			fits.writeto(dq_name,f_blob,overwrite=True)


main()








