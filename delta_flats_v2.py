#file: delta_flats_v2.py

import matplotlib.pyplot as plt
import csv
import pandas as pd 
import glob
from astropy.io import ascii 
import numpy as np 
import os
from astropy.io import fits
from scipy import ndimage
from photutils import CircularAnnulus
from photutils import CircularAperture



def read_in_bpix_098():
	filt_table = ascii.read('/user/holszewski/IR_flats/blob_sumry.txt', data_start=1, names= ['n','x','y','radius','flux','mjd','window','Flag'],delimiter=',' )
	#print(filt_table)blob_summary_new_copy.csv
	blobn = filt_table['n']
	blobx = filt_table['x']
	bloby = filt_table['y']
	blobrad = filt_table['radius']
	usaftermjd = filt_table['mjd']
	return (usaftermjd,blobx,bloby,blobrad,blobn)

def read_in_bpix():
	filt_table = ascii.read('/user/holszewski/IR_flats/blob_sumry_wide_edit_mjd.txt', data_start=1, names= ['n','x','y','radius','flux','mjd','app_win','flag'],delimiter=',')
	#print(filt_table)
	blobn = filt_table['n']
	blobx = filt_table['x']
	bloby = filt_table['y']
	blobrad = filt_table['radius']
	usaftermjd = filt_table['mjd']
	return (usaftermjd,blobx,bloby,blobrad,blobn)


def clac_deltas_error (im_name,mjd,blobn,blobx,bloby,blobrad,delta,day):
	sig = fits.getdata(im_name)
	sig[sig == 0] = 1
	im = 1/(sig**0.5)
	im[im == 1] = 0
	for i,date in enumerate(mjd):
		if date == day:
			num = blobn[i]
			#get the individual blob data at x-5,y-5
			x = float(blobx[i])-5
			y = float(bloby[i]) -5
			#print('1014',x,y)
			x_big = float(blobx[i])
			y_big = float(bloby[i])
			#print('1024',x_big,y_big)
			rad = (blobrad[i])
			x_min = int(x) -  (float(rad)+2)
			x_max = int(x) + (float(rad)+2)
			y_min = int(y) - (float(rad)+2)
			y_max = int(y) + (float(rad)+2)
			x_min_big = int(x_big) -  (float(rad)+2)
			x_max_big = int(x_big) + (float(rad)+2)
			y_min_big = int(y_big) - (float(rad)+2)
			y_max_big = int(y_big) + (float(rad)+2)
			if x_min < 0:
				x_min = 0
				x_min_big = 5
			if y_min < 0:
				y_min = 0
				y_min_big = 5
			if x_max > 1014:
				x_max = 1014
				x_max_big = 1019
			if y_max > 1014:
				y_max = 1014
				y_max_big = 1019

			aper = CircularAperture((x, y), rad)
			aper_big = CircularAperture((x_big, y_big), rad)

			I,J=np.meshgrid(np.arange(delta.shape[0]),np.arange(delta.shape[1]))
			#l,p=np.meshgrid(np.arange(im.shape[0]),np.arange(im.shape[1]))

			# calculate distance of all points to centre
			dist=np.sqrt((I-x_big)**2+(J-y_big)**2)
			#dist_small=np.sqrt((l-x)**2+(p-y)**2)
			
			
	
			sub_im = im[y_min:y_max,x_min:x_max]
			l,p=np.meshgrid(np.arange(sub_im.shape[0]),np.arange(sub_im.shape[1]))
			lx= rad+2
			ly=rad+2
			dist_small=np.sqrt((l-lx)**2+(p-ly)**2)
			#put the individual blob data into the 1024x1024 array at x,y
			#delta[y_min_big:y_max_big,x_min_big:x_max_big] = im[y_min:y_max,x_min:x_max]
			delta[np.where(dist<rad)]=sub_im[np.where(dist_small<rad)]
			result = ndimage.generic_filter(delta, np.nanmedian, size=3, mode='constant', cval=np.NaN)
			delta[delta==0] = result[delta==0]
			# Assign value of 1 to those points where dist<cr:



def main():
	list_filters = ['F160W','F098M','F105W','F125W','F140W','F110W']

	#read in table saying where blobs are and when
	for filt in list_filters:
		#if filt == 'F098M':
		#	mjd,blobx,bloby,blobrad,blobn = read_in_bpix_098()

		#else:
		mjd,blobx,bloby,blobrad,blobn = read_in_bpix()

		#create array of ones 1024x1024
		delta = np.zeros((1024,1024))
		
		current = os.getcwd()
		path = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_no512_run/' + filt +'/'#nov_4_no512_run
		os.chdir(path)
		#for each blob apearance date (a unique list of mjds) get the x,y of the blobs that apeard that day
		mjd_set = set(mjd)
		sorted_mjd=sorted(mjd_set)
		print(sorted_mjd)
		for day in sorted_mjd:
			#get the blobstaked image for that date
			#im_name = filt + '_' + str(day) + '_blob_median_test_intblobs.fits'#''_blob_median_v2.fits
			#error_name = filt + '_' + str(day) + '_blob_error_test_intblobs.fits'#''_blob_median_v2.fits
			sig_name = filt + '_' + str(day) + '_blob_signal_test_intblobs.fits'#''_blob_median_v2.fits
			clac_deltas_error(sig_name,mjd,blobn,blobx,bloby,blobrad,delta,day)
			#delta[delta==1] = 0
			#result = ndimage.generic_filter(delta, np.nanmean, size=3, mode='constant', cval=np.NaN)
			#delta[np.isnan(delta)] = 0



			output = filt + '_' + str(day) + '_delta_sig_error_ground_blobs_test3.fits'
			fits.writeto(output,delta,overwrite=True)

	os.chdir(path)
			
		
		
		
		
		
		
		
		
		
		
		
		
	
	#write out the 1024x1024 array as the detla for date
	
	
	
	
	

main()



















































	