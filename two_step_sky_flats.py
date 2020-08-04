#! /usr/bin/env python

"""Creats the IR sky flats.
Authors
-------
    Heather Kurtz, December 8, 2017
Use
-------

"""

import logging
import os
from multiprocessing import Pool
from bisect import bisect
from shutil import copyfile
import glob

import statistics
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clip
from astropy.io import ascii
from photutils import detect_sources
from photutils import detect_threshold

#from pyql.database.ql_database_interface import IR_flt_0
#from pyql.database.ql_database_interface import Master
#from pyql.database.ql_database_interface import session
from scipy import stats
from scipy.stats import sigmaclip
from wfc3tools import calwf3

from wfc3ir_tools import make_flattened_ramp_flt
from wfc3ir_tools import _reprocess_raw_crcorr
import earth_lim_cor 

# global variables: ie. logger
logger = logging.getLogger('testing')
hdlr = logging.FileHandler('/user/holszewski/IR_flats/noflat_mar5.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

filt_table = ascii.read('wide_filterlist.txt', data_start=1, delimiter=' ') #filter_info f110_info.txt 2nd_half_filters.txt wide_filterlist.txt
filt_list = filt_table['FILTER']
pflat_list = filt_table['TV3']
lflat_list = filt_table['Pipeline']
ff_list = [12065, 12068, 12453, 12460, 12451, 13495, 13496, 13498, 13504, 14037, 14038]
refilter = ['F105W','F110W']
asn_bad = ['idgb67010_asn.fits', 'idxm08010_asn.fits', 'idgbg7010_asn.fits', 'idgb93010_asn.fits',
            'idgb65010_asn.fits', 'idgb26010_asn.fits', 'idxm16010_asn.fits', 'idgbe8010_asn.fits', 
            'idxm13010_asn.fits', 'idgbk0010_asn.fits', 'idxm02010_asn.fits', 'idgb60010_asn.fits', 
            'idzf03010_asn.fits', 'idgb63010_asn.fits', 'idq203010_asn.fits', 
            'idgbf2010_asn.fits', 'idgbm5010_asn.fits', 'idgbl4010_asn.fits','idgb80010_asn.fits']
          #'ibohbb020_asn.fits','idix06010_asn.fits','ic8mb7030_asn.fits','icqtb7040_asn.fits','icqt93030_asn.fits','icqt94030_asn.fits','icqt82030_asn.fits',
          #'id8l03010_asn.fits','ic1602030_asn.fits','ic1602040_asn.fits','idia23020_asn.fits','idia23010_asn.fits','idia23030_asn.fits','ic8n40020_asn.fits',
          #'icqt60020_asn.fits','icwb0p020_asn.fits','idix02010_asn.fits','iboi6k020_asn.fits','iboi6o020_asn.fits','iboi6c020_asn.fits','iboi6f020_asn.fits',
          #'iboi6g020_asn.fits','iboi6h020_asn.fits','iboi6i020_asn.fits','iboi6j020_asn.fits','id8l03020_asn.fits','id5r10010_asn.fits','idix02020_asn.fits',
          #'idia22020_asn.fits','idia22030_asn.fits','idia22010_asn.fits','idia04010_asn.fits','idia04020_asn.fits','idia04030_asn.fits','ic8mb7040_asn.fits',
          #'icqtb5030_asn.fits','icqtb7030_asn.fits','icqt93040_asn.fits','icqt94040_asn.fits','icqt82040_asn.fits','ibohbe020_asn.fits', 'ibohbe020_asn.fits',
          #'idgb67010_asn.fits', 'idgb93010_asn.fits']


key_bad = ['/grp/hst/wfc3v/hkurtz/sky_flats/input_data/ibl71ki4q_raw.fits', '/grp/hst/wfc3v/hkurtz/sky_flats/input_data/iboi4qipq_raw.fits', 
           '/grp/hst/wfc3v/hkurtz/sky_flats/input_data/iboi4qihq_raw.fits', '/grp/hst/wfc3v/hkurtz/sky_flats/input_data/iboi4qieq_raw.fits',
           '/grp/hst/wfc3v/hkurtz/sky_flats/input_data/ibl71ki2q_raw.fits', '/grp/hst/wfc3v/hkurtz/sky_flats/input_data/iboi4qilq_raw.fits']
# TODO Create diagnostic plots for each FLT frame with a histogram of pixel values after masking. Overplot
# sigma-clipped mean and (sigma-clipped)
# TODO median. We can look at these later and decide which typically works better.
# TODO made dic to creat doc with /QL_path/file_flt.fits, filter, exptime, mean, median, mode, %Good_pix,
# (mean/med), (med-mode)/med,  Used?(y/n)


#def quary_ql(proid, filt):
#    # I need to use the ql quary to get: directories for images,
#
#    results = session.query(Master.dir, Master.rootname). \
#        join(IR_flt_0). \
#        filter(
#        # IR_flt_0.targname == 'tungsten',
#        IR_flt_0.filter == filt,
#        IR_flt_0.detector == 'ir',
#        # IR_flt_0.imagetyp == 'flat',
#        IR_flt_0.subarray == False,
#        IR_flt_0.exptime > 300,
#        IR_flt_0.proposid == proid)
#    # Turn the roots and dirs into locations we can use later.
#    locales = ['{}_flt.fits'.format(os.path.join(item.dir, item.rootname)) for item in results]
#
#    return locales


def open_file(file):
    hdulist = fits.open(file)
    data = hdulist[1].data
    #print(data)
    return data


def get_data(file):
    hdulist = fits.open(file)
    data = hdulist[1].data
    dq = hdulist[3].data
    hdulist.close()
    return data, dq


def open_update(file):
    hdulist = fits.open(file, mode = 'update')
    data = hdulist[1].data
    return data


def dq_convolved_mask(dq, data, ims):
    bit_mask = (4 + 8 + 32 + 128 + 512)
    dq0 = np.bitwise_and(dq, np.zeros(np.shape(dq), 'Int16') + bit_mask)
    dq0[dq0 > 0] = 1
    ims[dq0 > 0] = 1
    data[ims > 0] = np.nan


def dq_mask(dq, data):
    bit_mask = (4 + 8 + 32 + 128 + 512)
    dq0 = np.bitwise_and(dq, np.zeros(np.shape(dq), 'Int16') + bit_mask)
    dq0[dq0 > 0] = 1
    data[dq0 > 0] = np.mean(sigma_clip(data, sigma=5))


def weak_blob_mask(data):
    blob_mask = fits.getdata('/user/holszewski/IR_flats/test_ground_blobs.fits')
    blob_mask[blob_mask > 0] = 1
    data[blob_mask > 0] = np.mean(sigma_clip(data, sigma=5))


def weak_blob_mask_convolve(data):
    blob_mask = fits.getdata('/user/holszewski/IR_flats/test_ground_blobs.fits')
    blob_mask[blob_mask > 0] = 1
    data[blob_mask > 0] = np.nan


def find_sources(data):
    threshold = detect_threshold(data, snr=0.75)
    sigma = 1.5 * gaussian_fwhm_to_sigma  # FWHM = 2.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=4, filter_kernel=kernel)
    print(segm)
    seg = segm.data
    seg[seg > 0] = 1.0
    return seg


def persistince_source(file):
    try:
        p_file = file[:-9] + '_persist.fits'
        hdulist = fits.open(p_file)
        p_data = hdulist[1].data
    except FileNotFoundError:
        p_data = np.zeros((1014, 1014))
    return p_data


def persistince_masks(data, percist):
    data[percist > 0.005] = np.nan


def mask_sources(seg):
    seg[seg >= 0.012] = 1
    seg[seg < 0.012] = 0
    return seg


def convolve_data(seg_arr):
    g = Gaussian2DKernel(stddev=7)
    dataC = convolve(seg_arr, g, boundary='extend')
    return dataC


def write_file(file, hdr, plo, pro, end, filt):
    spro = str(pro)
    #file_name = '/grp/hst/wfc3v/hkurtz/sky_flats/short_run/' + filt + '/' + spro + '_' + file[-18:-8] + end
    file_name = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_run/' + filt + '/' + spro + '_' + file[-18:-8] + end
    prihdu = fits.PrimaryHDU(header=hdr)
    single_extension1 = fits.ImageHDU(data=plo.astype(np.float32))
    all_extensions = [prihdu, single_extension1]
    myhdulist = fits.HDUList(all_extensions)
    myhdulist.writeto(file_name, overwrite=True)


def normalize(data):
    values = data[~np.isnan(data)]
    values, clow, chigh = sigmaclip(values, low=3, high=3)
    mean = np.nanmean(values)
    image = data / mean
    return image, mean


def normalize_region(data):
    values = data[101:900, 101:900]
    nonan = values[~np.isnan(values)]
    value, clow, chigh = sigmaclip(nonan, low=5, high=5)
    mean = np.mean(value)
    image = data / mean
    return image, mean


#def stack(data_array_1):
#    image_median = np.nanmedian(data_array_1, axis=0)
#    image_mean = np.nanmean(data_array_1, axis=0)
#    image_sum = np.nansum(data_array_1, axis=0)
#    image_std = np.nanstd(data_array_1, axis=0)
#
#    return image_mean, image_median, image_sum, image_std


def sigclip(data):
    clip_data = sigma_clip(data, sigma=2, iters=3)
    return clip_data


def flat_field(filt, list_filt, tv3):
    for i in range(len(list_filt)):
        if list_filt[i] == filt:
            pflat = tv3[i]
            #lflat = pipeline[i]
        else:
            # print('flat_failed',filt)
            continue  # maybe add debugger here
    #n_file = '/grp/hst/cdbs/iref/' + pflat
    #o_file = '/grp/hst/cdbs/iref/' + lflat
    #o_flat = open_file(o_file)
    #n_flat = open_file(n_file)
    #o_data = o_flat[5:1019, 5:1019]
    #n_data = n_flat[5:1019, 5:1019]
    #new_data = data * (o_data / n_data)
    return pflat


def data_size(data):
    datsize = data.size
    nans = data[np.isnan(data)]
    nan_size = len(nans)
    return nan_size, datsize


def earth_lim_check(image):
    # nanim = image[~np.isnan(image)]
    # mode = statistics.mode(nanim)
    mean = np.nanmean(image)
    median = np.nanmedian(image)
    # diff = (median - mode)/median
    diff = mean / median
    return (diff)


def get_header_data(f):
    hdr1 = fits.getheader(f, 0)
    propid = hdr1['PROPOSID']
    calver = hdr1['CAL_VER']
    targ = hdr1['TARGNAME']
    expt = hdr1['EXPTIME']
    FILTER = hdr1['FILTER']

    logger.info('[PROPOSID:, CAL_VER:, TARGNAME:] %s', [propid, calver, targ])
    logger.info('EXPTIME: %s', expt)
    logger.info('FILTER: %s', FILTER)


def add_to_header(hdr, norm_mean, mean, median, std, good_pix, used):
    hdr['NORM'] = norm_mean
    hdr['Mean'] = mean
    hdr['Median'] = median
    hdr['STD'] = std
    hdr['PerGood'] = good_pix
    hdr['Used'] = used


def raw_header(file):
    hdr = fits.getheader(file, 0)
    samp = hdr['SAMP_SEQ']
    ap = hdr['APERTURE']
    filter_n = hdr['FILTER']
    date = hdr['EXPSTART']
    return(samp, ap, date,filter_n)


def read_bpix_file():
    filt_table = ascii.read('/user/holszewski/IR_flats/bpixtab_summary.txt', data_start=1, delimiter=',')
    bpixtab = filt_table['col1']
    usaftermjd = filt_table['col3']
    return (bpixtab,usaftermjd)


def match_date(date,useatermjd,bpixtab): 
    #dat = str(date)
    index = bisect(useatermjd, date) -1
    #corin=index -1
    match = bpixtab[index]
    match_name = match
    return(match_name)


def read_dark_file():
    filt_table = ascii.read('/user/holszewski/IR_flats/superdarks.txt', data_start=1, delimiter=',')
    dark = filt_table['superdark']
    sample = filt_table['ss']
    aper = filt_table['ap']
    mjd = filt_table['useafter']
    return (dark,sample,aper,mjd)


def read_missed_file():
    filt_table = ascii.read('/grp/hst/wfc3v/hkurtz/sky_flats/test/missed_file_lsit.txt', data_start=1, delimiter=' ')
    missed = filt_table['file']
    return (missed)


def match_dark(file,dark,samp,sample,dark_mjd,date):
    #dat = str(date)
    mjd_list = []
    index_list = []
    for i in range(len(sample)):
        #if (sample[i] == samp) and (aper[i] == ap):
        if sample[i] == samp.lower():
            mjd_list.append(dark_mjd[i])
            index_list.append(i)
        else:
            continue
    index = bisect(mjd_list, date) -1
    f_index = index_list[index]
    good_dark = dark[f_index]

    return(good_dark)


def segment_pipeline(f,propid,filter,ql_file):
    hdr1 = fits.getheader(f, 0)
    data, dq = get_data(f)
    print(f,'writing flt')
    write_file(f, hdr1, data, propid, 'flt.fits', filter)
    print(f, hdr1['EXPTIME'])
    #print('persistance for', ql_file)
    #p_data = persistince_source(ql_file)
    data_copy = np.copy(data)
    #data_mask[dq != 0] = 0
    dq_mask(dq, data_copy)
    weak_blob_mask(data_copy)
    #persistince_masks(data, percist)
    print('finding source')
    threshold = detect_threshold(data_copy, snr=0.75)
    sigma = 1.5 * gaussian_fwhm_to_sigma  # FWHM = 2.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    #kernel.normalize()
    segm = detect_sources(data_copy, threshold, npixels=4, filter_kernel=kernel)
    print(segm)
    seg = segm.data
    seg[seg > 0] = 1.0
    #seg = find_sources(data_mask)
    print(f,'writing seg')
    write_file(f, hdr1, seg, propid, 'seg.fits', filter)


def raw_2_flt(ql_file):
    current = os.getcwd()
    new_file = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + ql_file[-18:-8] + 'raw.fits'
    file_ql = ql_file[:-8] + 'raw.fits'
    copyfile(file_ql, new_file)
    print(new_file)
    hdulist = fits.open(new_file, mode = 'update')
    hdr = hdulist[0].header
    asn = hdr['ASN_TAB']
    if asn in asn_bad:
        print('no asn')
    else:
        asn_ql = ql_file[:-18] + asn
        asn_local = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + asn
        if asn == "NONE":
            print("have to make ASN")
        if asn != "NONE":
            try:
                copyfile(asn_ql, asn_local)
                print("ASN copied")
            except FileNotFoundError:
                print('no asn')
    samp,ap,date,filt = raw_header(new_file)
    bpix_value,useafter = read_bpix_file()
    match = match_date(date,useafter,bpix_value)
    hdr['BPIXTAB'] = match
    dark,sample,aper,dark_mjd = read_dark_file()
    g_dark = match_dark(new_file,dark,samp,sample,dark_mjd,date)
    hdr['DARKFILE'] = g_dark
    #tv3_flat = flat_field(filt, filt_list, pflat_list)
    #hdr['PFLTFILE'] = tv3_flat
    hdulist.close()
    flt = new_file[:-8] + 'flt.fits'
    ima = new_file[:-8] + 'ima.fits'
    if os.path.isfile(flt):
            os.remove(flt)
    if os.path.isfile(ima):
            os.remove(ima)
    #print("start calwf3")
    #os.chdir(new_file[:-18])
    #calwf3(new_file[-18:])
    #os.chdir(current)
    #print('finished calwf3 if there is no file we have a porblem')


def check_ff(ql_file, raw_f, lim = True):
    #raw_2_flt(ql_file)
    #current = os.getcwd()
    #os.chdir(raw_f[:-18])
    calwf3(raw_f[-18:])
    #os.chdir(current)
    print('finished calwf3 if there is no file we have a porblem')
    f=raw_f[:-8]+'flt.fits'
    print(f)
    logger.info(f)
    hdr1 = fits.getheader(f, 0)
    propid = hdr1['PROPOSID']
    calver = hdr1['CAL_VER']
    targ = hdr1['TARGNAME']
    expt = hdr1['EXPTIME']
    filter = hdr1['FILTER']
    asn = hdr1['ASN_TAB']
    

    logger.info('[PROPOSID:, CAL_VER:, TARGNAME:] %s', [propid, calver, targ])
    logger.info('EXPTIME: %s', expt)
    logger.info('FILTER: %s', filter)
    if lim == True:
        print('entering earth lim')
        earth_lim_cor.e_lim_cor(raw_f)
    if filter in refilter:
        #raw_file=f[:-8]+'raw.fits'
    #    print('entering cr_corr')
    #    print(raw_f)
        print('start asn')
        if asn == 'NONE': #make dummy asn table. remove this clause once calwf3 is patched
            new_asn_tab = np.rec.array([(os.path.basename(raw_f)[0:9],'EXP-DTH', 1)], formats = 'S14,S14,i1', names='MEMNAME,MEMTYPE,MEMPRSNT')
            hdu_1 = fits.BinTableHDU(new_asn_tab)
            none_hdu_list = fits.HDUList([fits.open(raw_f)[0], hdu_1])
            print('asn made')
            asn_name = raw_f[-18:-8] + 'asn.fits'
            hdu = fits.open(raw_f, mode='update')
            hdu[0].header['ASN_TAB'] = asn_name
            #hdr['ASN_TAB'] = asn_name
            hdu.close()
            print('write out asn')
            none_hdu_list.writeto(raw_f[:-8]+'asn.fits')

        elif asn in asn_bad: #make dummy asn table. remove this clause once calwf3 is patched
            new_asn_tab = np.rec.array([(os.path.basename(raw_f)[0:9],'EXP-DTH', 1)], formats = 'S14,S14,i1', names='MEMNAME,MEMTYPE,MEMPRSNT')
            hdu_1 = fits.BinTableHDU(new_asn_tab)
            none_hdu_list = fits.HDUList([fits.open(raw_f)[0], hdu_1])
            print('asn made')
            asn_name = raw_f[-18:-8] + 'asn.fits'
            hdu = fits.open(raw_f, mode='update')
            hdu[0].header['ASN_TAB'] = asn_name
            #hdr['ASN_TAB'] = asn_name
            hdu.close()
            print('write out asn')
            none_hdu_list.writeto(raw_f[:-8]+'asn.fits')

        else:
            asn_local = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + asn
            if os.path.isfile(asn_local):
                print('true')
        print('entering flattend ramp')
        make_flattened_ramp_flt(raw_f)#stats_subregion = None, outfile = None,
    #else:
        #continue
    #    print('filter_error')
#    if propid in ff_list:
#        if "PAR" in targ:
#            pipeline(f,propid,filter,ql_file)    
#        else:
#            print(f, "is a NOT a parallel observation")
#    elif propid not in ff_list:
#       pipeline(f,propid,filter,ql_file)
#    return


def pipeline(f,propid,filter,ql_file):

    logger.info(f)
    hdr1 = fits.getheader(f, 0)
    data, dq = get_data(f)
    print(f,'writing flt')
    write_file(f, hdr1, data, propid, 'flt.fits', filter)
    print(f, hdr1['EXPTIME'])
    print('persistance for', ql_file)
    p_data = persistince_source(ql_file)
    print('copy data')
    #data_copy = np.copy(data)
    ##data_mask[dq != 0] = 0
    #dq_mask(dq, data_copy)
    #weak_blob_mask(data_copy)
    #print('masked')
    ##read in segmentation map

    spro=str(propid)
    file_seg = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_run/' + filter + '/' + spro + '_' + f[-18:-8] + 'seg.fits'
    print(file_seg)
    seg = open_file(file_seg)
    print(seg)
    dataC = convolve_data(seg)
    write_file(f, hdr1, dataC, propid, 'con.fits', filter)
    im = mask_sources(dataC)
    write_file(f, hdr1, im, propid, 'cim.fits', filter)
    print('seg read in')
    dq_convolved_mask(dq, data, im)
    weak_blob_mask_convolve(data)
    write_file(f, hdr1, data, propid, 'dqm.fits', filter)
    print('DQ masked')
    print("per", p_data.shape, "data", data.shape)
    persistince_masks(data, p_data)
    print('persistance mask')
    write_file(f, hdr1, data, propid, 'per.fits', filter)
    image, norm_mean = normalize_region(data)
    print('file writed')
    write_file(f, hdr1, data, propid, 'norm.fits', filter)
    logger.info('Normalized to: %s', norm_mean)
    clipdata = sigclip(image)
    nan_siz, dat_siz = data_size(clipdata)
    #hdr1['NORM'] = norm_mean
    per = (nan_siz / dat_siz) * 100
    logger.info('Percent nan pixel: %s', (nan_siz / dat_siz))
    mean = np.nanmean(image)
    median = np.nanmedian(image)
    std = np.nanstd(image)
    print(f,norm_mean, mean, median, std, per)
    #write_file(f, hdr1, image, propid, 'tes.fits', filter)
    if nan_siz > (dat_siz * 0.75):
        logger.info('File has too many masked pixels. Not used.')
        used = 'per_bad'
        if ~np.isnan(mean) & ~np.isnan(median) & ~np.isnan(norm_mean) & ~np.isnan(std) & ~np.isnan(per):
            add_to_header(hdr1, norm_mean, mean, median, std, per, used)
        print(f,'writing pix')
        write_file(f, hdr1, image, propid, 'pix.fits', filter)


    # continue
    #else:
    #    diff = earth_lim_check(image)
    #    logger.info('Differance Mean-Median: %s', diff)
    #    if abs(diff) > 1.0011:
    #        logger.info('File has Earthlim. Not used.')
    #        used = 'Earthlim'
    #        if ~np.isnan(mean) & ~np.isnan(median) & ~np.isnan(norm_mean) & ~np.isnan(std) & ~np.isnan(per):
    #            add_to_header(hdr1, norm_mean, mean, median, std, per, used)
    #        print(f,'writing elm')
    #        write_file(f, hdr1, image, propid, 'elm.fits', filter)
    #    else:
    #        used = 'yes'
    #        if ~np.isnan(mean) & ~np.isnan(median) & ~np.isnan(norm_mean) & ~np.isnan(std) & ~np.isnan(per):
    #            add_to_header(hdr1, norm_mean, mean, median, std, per, used)
    #        print(f,'writing mdi')
    #        print(f, hdr1['EXPTIME'])
    #        write_file(f, hdr1, image, propid, 'mdi.fits', filter)
    #        #data_array[i, :, :] = clipdata
    #        logger.info('File Used.')
    else:
        used = 'yes'
        if ~np.isnan(mean) & ~np.isnan(median) & ~np.isnan(norm_mean) & ~np.isnan(std) & ~np.isnan(per):
            add_to_header(hdr1, norm_mean, mean, median, std, per, used)
        print(f,'writing mdi')
        print(f, hdr1['EXPTIME'])
        write_file(f, hdr1, image, propid, 'mdi.fits', filter)
        #data_array[i, :, :] = clipdata
        logger.info('File Used.')
    return ()


def run2(ql_file):
    raw_f = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + ql_file[-18:-8] + 'raw.fits'
    current = os.getcwd()
    os.chdir(raw_f[:-18])
    raw_2_flt(ql_file)
    check_ff(ql_file,raw_f)
    f=raw_f[:-8]+'flt.fits'
    hdr1 = fits.getheader(f, 0)
    propid = hdr1['PROPOSID']
    filter = hdr1['FILTER']
    targ = hdr1['TARGNAME']
    segment_pipeline(f,propid,filter,ql_file)
    #change flat
    hdu = fits.open(raw_f, mode='update')
    hdr = hdu[0].header
    hdr['FLATCORR'] = 'OMIT'
    hdu.close()
    #remove all flts and imas
    first_flts = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + ql_file[-18:-8] + 'flt.fits'
    first_imaflts = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/flattened_' + ql_file[-18:-8] + 'ima_flt.fits'
    first_imaimas = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/flattened_' + ql_file[-18:-8] + 'ima_ima.fits'    
    first_imas ='/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_input/' + ql_file[-18:-8] + 'ima.fits'
    existing_files = glob.glob(first_flts)
    existing_files.append(first_imas)
    existing_files.append(first_imaflts)
    existing_files.append(first_imaimas)

    print(existing_files)
    for file in existing_files:
        if os.path.isfile(file):
            os.remove(file)
            print('removed', file)

    print('running second time')
    check_ff(ql_file,raw_f,lim=False)
    if propid in ff_list:
        if "PAR" in targ:
            pipeline(f,propid,filter,ql_file)    
        else:
            print(f, "is a NOT a parallel observation")
    elif propid not in ff_list:
       pipeline(f,propid,filter,ql_file)
    os.chdir(current)







# move if statments as conditionals


# def sub_check(data):
#   ys,xs = np.shape(data)
#   if ys!=1014 or xs!=1014:
#       print ("Not a full array!")
#        return None,None


def main():

    # list_file=glob.glob('/grp/hst/wfc3a/GO_Links/12167/Visit05/*flt.fits')

    logger.info('The pipeline is starting')

    #pro_list = ['13000','14262','13667','14327','15118','11166','11101','11202','11208','11343','11557','11597','11600','11644',
    #            '11650','11666','11669','11838','11624','12051','11709','11738','12177','12194','12247','12265','12471','12487',
    #            '12496','12886','12942','12949','12990']

    #pro_list = ['11108', '11142', '11149', '11153', '11166', '11189', '11202', '11208', '11343', '11359',
    #            '11519', '11520', '11534', '11541', '11557', '11563', '11584', '11597', '11600', '11644',
    #            '11650', '11666', '11669', '11694', '11700', '11702', '11735', '11838', '11840', '11587',
    #            '11528', '11624', '12051', '12005', '11709', '11738', '11602', '11663', '12064', '13056',
    #            '12224', '12203', '11696', '12184', '12099', '12307', '12329', '12065', '12061', '12068',
    #            '12286', '12283', '12167', '12328', '12616', '12453', '12286', '12440', '12460',' 14038',
    #            '12581', '11636', '11734', '12060', '12062', '12063', '12064', '12177', '12194', '12247',
    #            '12265', '12442', '12443', '12444', '12445', '12471', '12487', '12496', '12498', '12451',
    #            '12578', '12764', '12886', '12905', '12930', '12942', '12949', '12959', '12960', '12974',
    #            '12990', '13000', '13002', '13045', '13110', '13117', '13294', '13303', '13480', '13614',
    #            '13641', '13644', '13688', '13718', '13792', '13793', '13831', '13844', '13868', '13951',
    #            '14262', '13667', '14327', '14459', '14699', '14718', '14719', '14721', '15118', '15137',
    #            '15287', '13495', '13496', '13498', '13504', '14037']
    list_files = []
   # for filt in filt_list:

   #     for i in range(len(pro_list)):
   #         logger.info('Getting the data for QL')
   #         list_file = quary_ql(pro_list[i], filt)
   #         for j in range(len(list_file)):
   #             list_files.append(list_file[j])
   #         print(len(list_files))
    current = os.getcwd()
    print('started')

    base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/MAST_data/'
    all_files = glob.glob(base_path+'*raw.fits')
    print(len(all_files))
    for f in all_files:
        #print(f)
        if f in key_bad:
            print(f,'bad_key')
        else:
            list_files.append(f)
        #    hdr = fits.getheader(f, 0)
        #    filt = hdr['FILTER']
        #    if filt == 'F105W':
        #        #new_file = '/grp/hst/wfc3v/hkurtz/sky_flats/el_105_input/' + f[-18:]
        #        list_files.append(f)
        #    elif filt == 'F140W':
        #        #new_file = '/grp/hst/wfc3v/hkurtz/sky_flats/el_105_input/' + f[-18:]
        #        list_files.append(f)
        #    else:
        #        continue
    print(len(list_files))

    p = Pool(30)
    result = p.map(run2, list_files)
    print(result)


main()
