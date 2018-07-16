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
import multiprocessing
from multiprocessing import Pool
import os
import sys
from typing import Tuple

from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import convolve
from astropy.stats import sigma_clip
from numpy.core.multiarray import ndarray
from scipy import stats
# from mx import DateTime
import numpy as np
from photutils import DAOStarFinder
from photutils import detect_sources
from photutils import detect_threshold
# import Statistics
# import Stats2
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


def quary_ql(proid, filt):
    # I need to use the ql quary to get: directories for images,

    results = session.query(Master.dir, Master.rootname). \
        join(IR_flt_0). \
        filter(
        # IR_flt_0.targname == 'tungsten',
        IR_flt_0.filter == filt,
        IR_flt_0.detector == 'ir',
        # IR_flt_0.imagetyp == 'flat',
        IR_flt_0.subarray == False,
        IR_flt_0.exptime > 600,
        IR_flt_0.proposid == proid)
    # Turn the roots and dirs into locations we can use later.
    locales = ['{}_flt.fits'.format(os.path.join(item.dir, item.rootname)) for item in results]

    return locales


def get_data(file):
    hdulist = fits.open(file)
    data = hdulist[1].data
    # masked_data=hdulist[1].data
    dq = hdulist[3].data
    hdulist.close()
    return (data, dq)


def dq_mask(dq, data, ims):
    bit_mask = (4 + 16 + 32 + 128 + 512)
    dq0 = np.bitwise_and(dq, np.zeros(np.shape(dq), 'Int16') + bit_mask)
    dq0 == 0
    dq0[dq0 > 0] = 1
    # data[dq!=0]=np.nan
    ims[dq0 > 0] = 1
    data[ims > 0] = np.nan


def find_sources(data):
    threshold = detect_threshold(data, snr=0.85)
    sigma = 2.0 * gaussian_fwhm_to_sigma  # FWHM = 2.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    print(kernel)
    segm = detect_sources(data, threshold, npixels=10, filter_kernel=kernel)
    return (segm)


def persistince_masks(data, per):
    data[per > 0.005] = np.nan
    return data


def mask_sources(seg):
    seg[seg >= 0.010] = 1
    seg[seg < 0.010] = 0
    return (seg)


def convolv_data(seg_arr):
    g = Gaussian2DKernel(stddev=7)
    # Convolve data
    dataC = convolve(seg_arr, g, boundary='extend')
    return (dataC)


def write_file(file, hdr, plo, pro):
    spro = str(pro)
    file_name = '/user/hkurtz/IR_flats/test_f098/' + spro + file[-18:-8] + 'mdi.fits'
    prihdu = fits.PrimaryHDU(header=hdr)
    single_extension1 = fits.ImageHDU(data=plo.astype(np.float32))
    all_extensions = [prihdu, single_extension1]
    myhdulist = fits.HDUList(all_extensions)
    myhdulist.writeto(file_name, overwrite=True)


def write_seg(file, hdr, plo, pro):
    spro = str(pro)
    file_name = '/user/hkurtz/IR_flats/test_f098/' + spro + file[-18:-8] + 'seg.fits'
    prihdu = fits.PrimaryHDU(header=hdr)
    single_extension1 = fits.ImageHDU(data=plo.astype(np.float32))
    all_extensions = [prihdu, single_extension1]
    myhdulist = fits.HDUList(all_extensions)
    myhdulist.writeto(file_name, overwrite=True)


def normalize(data):  # use region [101:900,101:900]
    values = data[~np.isnan(data)]
    values, clow, chigh = sigmaclip(values, low=3, high=3)
    mean = np.mean(values)
    image = data / mean
    return (image, mean)


def norm_region(data):  # use region [101:900,101:900]
    values = data[101:900, 101:900]
    values, clow, chigh = sigmaclip(values, low=3, high=3)
    mean = np.mean(values)
    image = data / mean
    return (image, mean)


def Stack(data_array_1):
    image_median = np.nanmedian(data_array_1, axis=0)
    image_mean = np.nanmean(data_array_1, axis=0)
    image_sum = np.nansum(data_array_1, axis=0)
    image_std = np.nanstd(data_array_1, axis=0)

    return (image_mean, image_median, image_sum, image_std)


def sigclip(data):
    # sigmaclip the data
    clip_data = sigma_clip(data, sigma=2, iters=3)
    return (clip_data)


def data_size(data):
    datsize = data.size
    nans = data[np.isnan(data)]
    nan_size = len(nans)
    return (nan_size, datsize)


def earth_lim_check(image):
    mode = stats.mode(~np.isnan(image))
    #mean = np.nanmean(image)
    median = np.nanmedian(image)
    diff = (median - mode)/median


def testing(f):
    logger = multiprocessing.get_logger()
    hdlr = logging.FileHandler('/user/hkurtz/IR_flats/098.log')
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    logger.info(f)
    print(f)
    hdr1 = fits.getheader(f, 0)
    propid = hdr1['PROPOSID']
    calver = hdr1['CAL_VER']
    targ = hdr1['TARGNAME']
    expt = hdr1['EXPTIME']
    FILTER = hdr1['FILTER']

    logger.info('[PROPOSID:, CAL_VER:, TARGNAME:] %s', [propid, calver, targ])
    # logger.info('CAL_VER: %s',calver)
    # logger.info('TARGNAME: %s',targ)
    logger.info('EXPTIME: %s', expt)
    logger.info('FILTER: %s', FILTER)
    data, dq = get_data(f)
    data_mask = np.copy(data)
    data_mask[dq != 0] = 0
    # fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    # ax1.imshow(data_mask, origin='lower')
    # print('data',data[200:210,200:210])
    # print('data masked',data_mask[200:210,200:210])
    segm = find_sources(data_mask)
    # fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
    # ax1.imshow(segm, origin='lower',cmap=segm.cmap(random_state=12345))
    # plt.show()

    seg = segm.array  # add to find source function for the next 2 lines
    seg[seg > 0] = 1.0
    write_seg(f, hdr1, seg, propid)
    dataC = convolv_data(seg)
    im = mask_sources(dataC)
    dq_mask(dq, data, im)
    image, norm_mean = normalize(data)
    logger.info('Normalized to: %s', norm_mean)
    clipdata = sigclip(image)
    nan_siz, dat_siz = data_size(clipdata)
    hdr1['NORM'] = norm_mean
    logger.info('Percent nan pixel: %s', (nan_siz / dat_siz))
    if nan_siz > (dat_siz * 0.8):
        # list_bad.append(f)
        logger.info('File has too many masked pixels. Not used.')
    # continue
    else:
        mean = np.nanmean(image)
        median = np.nanmedian(image)
        diff = mean - median
        logger.info('Differance Mean-Median: %s', diff)
        if abs(diff) > 1.0:
            # list_lim.append(f)
            logger.info('File has Earthlim. Not used.')
        else:
            write_file(f, hdr1, image, propid)
            # data_array[i, :, :] = clipdata
            # list_good.append(f)
            logger.info('File Used.')


# move if statments as conditionals


# def sub_check(data):
#	ys,xs = np.shape(data)           
#	if ys!=1014 or xs!=1014:
#		print ("Not a full array!")
#        return None,None


def main():
    # quarry ql put values into dictionary
    # put into forloop for filters

    # for i in len(dic[rootnames]):
    #	direct=dic[dir(i)]

    # open file

    # test file to see it good to use

    #	if quality =='good':
    # do everything else (find source, mask source normalize, and save)
    #	if qality=='bad':
    #		continue

    # test data:

    # list_file=glob.glob('/grp/hst/wfc3a/GO_Links/12167/Visit05/*flt.fits')

    logger = logging.getLogger('f098')
    hdlr = logging.FileHandler('/user/hkurtz/IR_flats/f098.log')
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    # logging.basicConfig(level=logging.INFO)
    # logger = logging.getLogger(__name__)
    logger.info('The pipeline is starting')

    # pro_list = ['11702']
    pro_list = ['11108', '11142', '11149', '11153', '11166', '11189', '11202', '11208', '11343', '11359',
                '11519', '11520', '11534', '11541', '11557', '11563', '11584', '11597', '11600', '11644',
                '11650', '11666', '11669', '11694', '11700', '11702', '11735', '11838', '11840', '11587']
    list_files = []
    for i in range(len(pro_list)):
        logger.info('Getting the data for QL')
        list_file = quary_ql(pro_list[i], 'F098M')
        for j in range(len(list_file)):
            list_files.append(list_file[j])
        # print(len(list_files))
    # list_file=quary_ql('12025','F160W')
    hdr = fits.getheader(list_files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf = len(list_files)
    set_data = fits.getdata(list_files[0], 1)
    # data_array = np.empty((nf, ny, nx), dtype=float)
    list_bad = []
    list_good = []
    list_lim = []
    p = Pool(8)
    result = p.map(testing, list_files)
    print(result)


# for i , f in enumerate(list_files):
#		logger.info(f)
#		print(f)
#		hdr1 = fits.getheader(f, 0)
#		propid=hdr1['PROPOSID']
#		calver=hdr1['CAL_VER']
#		targ=hdr1['TARGNAME']
#		expt=hdr1['EXPTIME']
#		FILTER=hdr1['FILTER']
#		logger.info('PROPOSID: %s', propid)
#		logger.info('CAL_VER: %s',calver)
#		logger.info('TARGNAME: %s',targ)
#		logger.info('EXPTIME: %s',expt)
#		logger.info('FILTER: %s',FILTER)
#		data,dq=get_data(f)
#		data_mask=np.copy(data)
#		data_mask[dq!=0]=0
#		#fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
#		#ax1.imshow(data_mask, origin='lower')
#		#print('data',data[200:210,200:210])
#		#print('data masked',data_mask[200:210,200:210])
#		segm=find_sources(data_mask)
#		#fig, (ax1) = plt.subplots(1, 1, figsize=(5, 5))
#		#ax1.imshow(segm, origin='lower',cmap=segm.cmap(random_state=12345))
#		#plt.show()
#
#		seg=segm.array
#		seg[seg>0]=1.0
#		dataC=convolv_data(seg)
#		im=mask_sources(dataC)
#		dq_mask(dq,data,im)
#		image,norm_mean=normalize(data)
#		logger.info('Normalized to: %s', norm_mean)
#		clipdata=sigclip(image)
#		nan_siz,dat_siz=data_size(clipdata)
#		hdr1['NORM']=norm_mean
#		logger.info('Number of nan pixel: %s', nan_siz)
#		if nan_siz>(dat_siz*0.8):
#			list_bad.append(f)
#			logger.info('File has too many masked pixels. Not used.')
#			#continue
#		else:
#			mean=np.nanmean(image)
#			median=np.nanmedian(image)
#			diff=mean-median
#			logger.info('Differance Mean-Median: %s',diff)
#			if abs(diff)>1.0:
#				list_lim.append(f)
#				logger.info('File has Earthlim. Not used.')
#			else:
#				#wrtie_file(f,hdr1,image,propid)
#				#data_array[i, :, :] = clipdata
#				list_good.append(f)
#				logger.info('File Used.')
# print('good',len(list_good))
# print('bad',len(list_bad))
# print('lim',len(list_lim))
# S_mean,S_median,S_sum,S_std=Stack(data_array)

# fits.writeto('test_F160_10_mean_sn5.fits', S_mean,overwrite=True)
# fits.writeto('test_F160_10_median_sn5.fits', S_median,overwrite=True)


main()
