#! /usr/bin/env python

"""Creats the IR sky flats.
Authors
-------
    Heather Kurtz, December 8, 2017
Use
-------

"""

import glob
import logging
from multiprocessing import Pool
import os
import sys

from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import convolve
#from mx import DateTime
import numpy as np
from photutils import DAOStarFinder
from photutils import detect_sources
from photutils import detect_threshold
#import Statistics
#import Stats2
import tempfile
import threading
from threading import Thread
import scipy
from scipy.stats import sigmaclip
from scipy.ndimage import gaussian_filter
from scipy import ndimage



from pyql.database.ql_database_interface import session
from pyql.database.ql_database_interface import Master
from pyql.database.ql_database_interface import IR_flt_0


def quary_ql():
#I need to use the ql quary to get: directories for images,


    results = session.query(Master.dir, Master.rootname).\
          join(IR_flt_0).\
          filter(
              #IR_flt_0.targname == 'tungsten',
              #IR_flt_0.filter == 'f105w',
              IR_flt_0.detector == 'ir',
              #IR_flt_0.imagetyp == 'flat',
              #IR_flt_0.exptime > 1,
              IR_flt_0.proposid == '11528')
# Turn the roots and dirs into locations we can use later.
    locales = ['{}_flt.fits'.format(os.path.join(item.dir, item.rootname)) for item in results]

    return locales

def get_data(file):

    hdulist=fits.open(file)
    data=hdulist[1].data
    dq=hdulist[3].data
    hdulist.close()
    return( data, dq)


def dq_mask(dq,data,ims):
    bit_mask = (4+16+32+128+512)
    dq0 = np.bitwise_and(dq,np.zeros(np.shape(dq),'Int16')+ bit_mask)
    dq0==0
    dq0[dq0>0]=1
    data[dq!=0]=np.nan
    ims[dq0>0]=1
    data[ims>0]=np.nan


def find_sources(data):

    threshold = detect_threshold(data, snr=3.)
    sigma = 2.0 * gaussian_fwhm_to_sigma    # FWHM = 2.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    print(kernel)
    segm = detect_sources(data, threshold, npixels=5)#, filter_kernel=kernel)
    return(segm)
    

#def persistince_masks():



def mask_sources(seg):	
  maps=seg.array
  maps[maps>0]=1000
  sigma=10.0 * gaussian_fwhm_to_sigma
  im=scipy.ndimage.gaussian_filter(maps, sigma, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)
  im[im>0]=1
  return(im)




def convolv_data(data):
	k=np.ones((3,3))
	dataC=convolve(data, k, boundary='fill', fill_value=1.0, nan_treatment='interpolate')
	return(dataC)

# clip,low,high=sigmaclip(data)
# s=clip.std()

# k = np.array([[1,1,1],[1,1,1],[1,1,1]])
# dataC=ndimage.convolve(data, k, mode='constant', cval=0.0)


def wrtie_file(file,plo):
  file_name=file[38:-8]+'mdi.fits'
  fits.writeto(file_name, plo,overwrite=True)


def normalize(data):
  values = data[~np.isnan(data)]
  values, clow, chigh = sigmaclip(values, low=3, high=3)
  mean = np.mean(values)
  image=data/mean
  return(image)


def Stack(data_array_1):
 #hdr = fits.getheader(list_of_files[0], 1)
 #nx = hdr['NAXIS1']
 #ny = hdr['NAXIS2']
 #nf = len(list_of_files)
 #set_data=fits.getdata(list_of_files[0], 1)
 ##makes empty array to be filled with data
 #data_array = np.empty((nf, ny, nx), dtype=float)
 #for i , f in enumerate(list_of_files):
 #  data_1=fits.getdata(f, 1)
 #  data_array[i, :, :] = data_1
  image_median_1 = np.median(data_array_1, axis=0)
  image_mean_1= np.mean(data_array_1, axis=0)
  return(image_mean_1,image_median_1)




def main():
	#quarry ql put values into dictionary
	#put into forloop for filters

	#for i in len(dic[rootnames]):
	#	direct=dic[dir(i)]
				
		#open file

		#test file to see it good to use

	#	if quality =='good':
			#do everything else (find source, mask source normalize, and save)
	#	if qality=='bad':
	#		continue


	#test data:
  
  list_file=glob.glob('/grp/hst/wfc3a/GO_Links/12167/Visit05/*flt.fits')
  hdr = fits.getheader(list_file[0], 1)
  nx = hdr['NAXIS1']
  ny = hdr['NAXIS2']
  nf = len(list_file)
  set_data=fits.getdata(list_file[0], 1)
  data_array = np.empty((nf, ny, nx), dtype=float)
  for i , f in enumerate(list_file):
    data,dq=get_data(f)
    #find_sources(data)
    dataC=convolv_data(data)
    segm=find_sources(dataC)
    im=mask_sources(segm)
    dq_mask(dq,data,im)
    image=normalize(data)
    wrtie_file(f,data)
    data_array[i, :, :] = image
  S_mean,S_median=Stack(data_array)
  fits.writeto('stacked_mean.fits', S_mean,overwrite=True)
  fits.writeto('stacked_median.fits', S_median,overwrite=True)



main()

  

	#read in masked images for stacking
	#combine all images






















