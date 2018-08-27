#! /usr/bin/env python

"""Creats the IR sky flats.
Authors
-------
    Heather Kurtz, December 8, 2017
Use
-------

"""

import logging
import multiprocessing
import os
from multiprocessing import Pool

import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clip
from astropy.io import ascii
from photutils import detect_sources
from photutils import detect_threshold
from pyql.database.ql_database_interface import IR_flt_0
from pyql.database.ql_database_interface import Master
from pyql.database.ql_database_interface import session
from scipy import stats
from scipy.stats import sigmaclip

# global variables: ie. logger
logger = logging.getLogger('test')
hdlr = logging.FileHandler('/user/hkurtz/IR_flats/test.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)                                                 
logger.setLevel(logging.INFO)

filt_table = ascii.read('filter_info', data_start=1 , delimiter=' ')
filt_list = filt_table['FILTER']
pflat_list = filt_table['TV3']
lflat_list = filt_table['Pipeline']

# TODO Create diagnostic plots for each FLT frame with a histogram of pixel values after masking. Overplot  sigma-clipped mean and (sigma-clipped)
# TODO median. We can look at these later and decide which typically works better.
# TODO made dic to creat doc with /QL_path/file_flt.fits, filter, exptime, mean, median, mode, %Good_pix, (mean/med), (med-mode)/med,  Used?(y/n)


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


def open_file(file):
    hdulist = fits.open(file)
    data = hdulist[1].data
    return data


def get_data(file):
    data = open_file(file)
    # masked_data=hdulist[1].data
    dq = hdulist[3].data
    hdulist.close()
    return data, dq


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
    seg = segm.array
    seg[seg > 0] = 1.0
    return seg


def persistince_source(file):
    p_file=file[:-19]+'Persist/'+file[-19:-9]+'persist.fits'
    hdulist = fits.open(p_file)
    p_data = hdulist[1].data
    return p_data


def persistince_masks(data, percist):
    data[percist > 0.005] = np.nan


def mask_sources(seg):
    seg[seg >= 0.010] = 1
    seg[seg < 0.010] = 0
    return seg


def convolve_data(seg_arr):
    g = Gaussian2DKernel(stddev=7)
    # Convolve data
    dataC = convolve(seg_arr, g, boundary='extend')
    return dataC


def write_file(file, hdr, plo, pro, end):
    spro = str(pro)
    file_name = '/user/hkurtz/IR_flats/test_f098/' + spro + file[-18:-8] + end
    prihdu = fits.PrimaryHDU(header=hdr)
    single_extension1 = fits.ImageHDU(data=plo.astype(np.float32))
    all_extensions = [prihdu, single_extension1]
    myhdulist = fits.HDUList(all_extensions)
    myhdulist.writeto(file_name, overwrite=True)


def normalize(data):
    values = data[~np.isnan(data)]
    values, clow, chigh = sigmaclip(values, low=3, high=3)
    mean = np.mean(values)
    image = data / mean
    return image, mean


def normalize_region(data):  # use region [101:900,101:900]
    values = data[101:900, 101:900]
    values, clow, chigh = sigmaclip(values, low=3, high=3)
    mean = np.mean(values)
    image = data / mean
    return image, mean


def stack(data_array_1):
    image_median = np.nanmedian(data_array_1, axis=0)
    image_mean = np.nanmean(data_array_1, axis=0)
    image_sum = np.nansum(data_array_1, axis=0)
    image_std = np.nanstd(data_array_1, axis=0)

    return image_mean, image_median, image_sum, image_std


def sigclip(data):
    # sigmaclip the data
    clip_data = sigma_clip(data, sigma=2, iters=3)
    return clip_data

def flat_field(data, filter, list_filt, tv3, pipeline):
    for i in list_filt:
        if list_filt[i] == filter:
            pflat = tv3[i]
            lflat = pipeline[i]
        else:
            continue
    n_file = '/grp/hst/cdbs/iref/' + pflat
    o_file = '/grp/hst/cdbs/iref/' + lflat
    o_flat = open_file(o_file)
    n_flat = open_file(n_file)
    o_data=o_flat[6:1019,6:1019]
    n_data=n_flat[6:1019,6:1019]
    new_data=data*(o_data/n_data)
    return new_data


def data_size(data):
    datsize = data.size
    nans = data[np.isnan(data)]
    nan_size = len(nans)
    return nan_size, datsize


def earth_lim_check(image):
    mode = stats.mode(~np.isnan(image))
    #mean = np.nanmean(image)
    median = np.nanmedian(image)
    diff = (median - mode)/median
    return(diff,mode)


def get_header_data(f):
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


def add_to_header(hdr,norm_mean, mean, median, mode, good_pix, used):
    hdr['NORM'] = norm_mean
    hdr['Mean'] = mean
    hdr['Median'] = median
    hdr['Mode'] = mode
    hdr['PerGood'] = good_pix
    hdr['Used'] = used


def testing(f):
    #logger = multiprocessing.get_logger()
    #hdlr = logging.FileHandler('/user/hkurtz/IR_flats/098.log')
    #formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    #hdlr.setFormatter(formatter)
    #logger.addHandler(hdlr)
    #logger.setLevel(logging.INFO)

    logger.info(f)
    print(f)
    hdr1 = fits.getheader(f, 0)
    propid = hdr1['PROPOSID']
    calver = hdr1['CAL_VER']
    targ = hdr1['TARGNAME']
    expt = hdr1['EXPTIME']
    filter = hdr1['FILTER']

    logger.info('[PROPOSID:, CAL_VER:, TARGNAME:] %s', [propid, calver, targ])
    # logger.info('CAL_VER: %s',calver)
    # logger.info('TARGNAME: %s',targ)
    logger.info('EXPTIME: %s', expt)
    logger.info('FILTER: %s', filter)
    data_pipe, dq = get_data(f)
    data = flat_field(data_pipe, filter, filt_list, pflat_list, lflat_list)
    p_data = persistince_source(f)
    data_mask = np.copy(data)
    data_mask[dq != 0] = 0
    seg = find_sources(data_mask)
    write_file(f, hdr1, seg, propid, 'seg.fits')
    dataC = convolve_data(seg)
    im = mask_sources(dataC)
    dq_mask(dq, data, im)
    persistince_masks(data, p_data)
    image, norm_mean = normalize_region(data)
    logger.info('Normalized to: %s', norm_mean)
    clipdata = sigclip(image)
    nan_siz, dat_siz = data_size(clipdata)
    hdr1['NORM'] = norm_mean
    per = (nan_siz / dat_siz)*100
    logger.info('Percent nan pixel: %s', (nan_siz / dat_siz))
    if nan_siz > (dat_siz * 0.75):
        # list_bad.append(f)
        logger.info('File has too many masked pixels. Not used.')
    # continue
    else:
        mean = np.nanmean(image)
        median = np.nanmedian(image)
        diff,mode = earth_lim_check(image)
        logger.info('Differance Mean-Median: %s', diff)
        if abs(diff) > 1.0011:
            # list_lim.append(f)
            logger.info('File has Earthlim. Not used.')
        else:
            used = 'yes'

            add_to_header(hdr1,norm_mean, mean, median, mode, per, used)
            write_file(f, hdr1, image, propid, 'mdi.fits')
            # data_array[i, :, :] = clipdata
            # list_good.append(f)
            logger.info('File Used.')
    return()

# move if statments as conditionals


# def sub_check(data):
#	ys,xs = np.shape(data)           
#	if ys!=1014 or xs!=1014:
#		print ("Not a full array!")
#        return None,None


def main():
    # list_file=glob.glob('/grp/hst/wfc3a/GO_Links/12167/Visit05/*flt.fits')

    logger.info('The pipeline is starting')



    pro_list = ['11702']
    #pro_list = ['11108', '11142', '11149', '11153', '11166', '11189', '11202', '11208', '11343', '11359',
     #           '11519', '11520', '11534', '11541', '11557', '11563', '11584', '11597', '11600', '11644',
      #          '11650', '11666', '11669', '11694', '11700', '11702', '11735', '11838', '11840', '11587',
    #            '11528', '11624', '12051', '12005', '11709', '11738', '11602', '11663', '12064', '12197',
    #            '12224', '12203', '11696', '12184', '12099', '12307', '12329', '12065', '12061', '12068',
    #            '12286', '12283', '12167', '11591', '12328', '12616', '12453', '12286', '12440', '12460',
    #            '12581', '11636', '11734', '12060', '12062', '12063', '12064', '12177', '12194', '12247',
    #            '12265', '12442', '12443', '12444', '12445', '12471', '12487', '12496', '12498', '12451',
    #            '12578', '12764', '12886', '12905', '12930', '12942', '12949', '12959', '12960', '12974',
    #            '12990', '13000', '13002', '13045', '13110', '13117', '13294', '13303', '13480', '13614',
    #            '13641', '13644', '13688', '13718', '13792', '13793', '13831', '13844', '13868', '13951',
    #            '14262', '13667', '14327', '14459', '14699', '14718', '14719', '14721', '15118', '15137',
    #            '15287', '13495', '13496', '13498', '13504', '14307', '14308']
    list_files = []
    for filt in filt_list['filter']:

        for i in range(len(pro_list)):
            logger.info('Getting the data for QL')
            list_file = quary_ql(pro_list[i], filt)
            for j in range(len(list_file)):
                list_files.append(list_file[j])
            # print(len(list_files))
        # list_file=quary_ql('12025','F160W')
        #hdr = fits.getheader(list_files[0], 1)
        #nx = hdr['NAXIS1']
        #ny = hdr['NAXIS2']
        #nf = len(list_files)
        #set_data = fits.getdata(list_files[0], 1)
        # data_array = np.empty((nf, ny, nx), dtype=float)
        p = Pool(8)
        result = p.map(testing, list_files)
        print(result)


# print('good',len(list_good))
# print('bad',len(list_bad))
# print('lim',len(list_lim))
# S_mean,S_median,S_sum,S_std=Stack(data_array)

# fits.writeto('test_F160_10_mean_sn5.fits', S_mean,overwrite=True)
# fits.writeto('test_F160_10_median_sn5.fits', S_median,overwrite=True)


main()
