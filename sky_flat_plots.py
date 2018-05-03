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
from astropy.stats import sigma_clip

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
import matplotlib.pyplot as plt



from pyql.database.ql_database_interface import session
from pyql.database.ql_database_interface import Master
from pyql.database.ql_database_interface import IR_flt_0


def quary_ql(proid,filt):
#I need to use the ql quary to get: directories for images,


    results = session.query(Master.dir, Master.rootname).\
          join(IR_flt_0).\
          filter(
              #IR_flt_0.targname == 'tungsten',
              IR_flt_0.filter == filt,
              IR_flt_0.detector == 'ir',
              #IR_flt_0.imagetyp == 'flat',
              IR_flt_0.exptime > 300,
              IR_flt_0.proposid == proid)
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
    bit_mask = (4+16+32+128)#+512)
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
    segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
    return(segm)
    

#def persistince_masks():


def mask_sources(seg):
	seg[seg>0.003]=1
	seg[seg<0.003]=0
	return(seg)


def convolv_data(seg_arr):
	g = Gaussian2DKernel(stddev=15)
	# Convolve data
	dataC = convolve(seg_arr, g, boundary='extend')
	return(dataC)


def wrtie_file(file,plo,pro):
  spro=str(pro)
  file_name=spro+file[-18:-8]+'mdi.fits'
  fits.writeto(file_name, plo,overwrite=True)


def normalize(data):
  values = data[~np.isnan(data)]
  values, clow, chigh = sigmaclip(values, low=3, high=3)
  mean = np.mean(values)
  image=data/mean
  return(image)


def Stack(data_array_1):
  image_median = np.nanmedian(data_array_1, axis=0)
  image_mean= np.nanmean(data_array_1, axis=0)
  image_sum= np.nansum(data_array_1, axis=0)
  image_std= np.nanstd(data_array_1, axis=0)

  return(image_mean,image_median,image_sum,image_std)


def sigclip(data):
  #sigmaclip the data
  clip_data = sigma_clip(data, sigma=2, iters=3)
  return(clip_data)

def data_size(data):
  nans = np.isnan(data)
  no_nans_data = data[~nans]
  nan_size=no_nans_data.size
  datsize=data.size
  return(nan_size,datsize)
  

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
  
  #list_file=glob.glob('/grp/hst/wfc3a/GO_Links/12167/Visit05/*flt.fits')
  pro_list=['12167']
  #pro_list=['11108','11142','11149','11153','11166','11189','11202','11208','11343','11359',
   #         '11519','11520','11534','11541','11557','11563','11584','11597','11600','11644',
    #        '11650','11666','11669','11694','11700','11702','11735','11838','11840','11587']
  list_files=[]
  for i in range(len(pro_list)):
    list_file=quary_ql(pro_list[i],'F160W')
    for j in range(len(list_file)):
      list_files.append(list_file[j])
  #list_file=quary_ql('12025','F160W')
  list_files=['/grp/hst/wfc3a/GO_Links/12167/Visit02/ibhg02q2q_flt.fits','/grp/hst/wfc3a/GO_Links/12167/Visit10/ibhg10b6q_flt.fits']
  hdr = fits.getheader(list_files[0], 1)
  nx = hdr['NAXIS1']
  ny = hdr['NAXIS2']
  nf = len(list_files)
  set_data=fits.getdata(list_files[0], 1)
  data_array = np.empty((nf, ny, nx), dtype=float)
  for i , f in enumerate(list_files):
    print(f)
    hdr1 = fits.getheader(f, 0)
    propid=hdr1['PROPOSID']
    data,dq=get_data(f)
    segm=find_sources(data)
    fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    ax1.imshow(segm, origin='lower',cmap=segm.cmap(random_state=12345))
    plt.show()
    seg=segm.array
    seg[seg>0]=1.0
    dataC=convolv_data(seg)
    im=mask_sources(dataC)
    fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    ax1.imshow(im, origin='lower')
    plt.show()
    dq_mask(dq,data,im)
    image=normalize(data)
    fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    ax1.imshow(image, origin='lower')
    plt.show()
    clipdata=sigclip(image)
    fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    ax1.imshow(clipdata, origin='lower')
    plt.show()
    nan_siz,dat_siz=data_size(clipdata)
    wrtie_file(f,image,propid)
    if nan_siz>(dat_siz*0.8):
        continue
    else:
        wrtie_file(f,image,propid)
        data_array[i, :, :] = clipdata
  #S_mean,S_median,S_sum,S_std=Stack(data_array)

  #fits.writeto('test_F160_10_mean_sn5.fits', S_mean,overwrite=True)
  #fits.writeto('test_F160_10_median_sn5.fits', S_median,overwrite=True)



main()

  

	#read in masked images for stacking
	#combine all images






















