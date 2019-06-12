#! /usr/bin/env python
#file: earth_lim_cor.py

import matplotlib.pyplot as plt
import pstat
from wfc3tools import calwf3
import numpy as np
import os
from astropy.io import fits


def e_lim_cor(file):
	print(file)
	ima = file[:-8] + 'ima.fits'
	flt = file[:-8] + 'flt.fits'
	ratio = get_ratio(ima)
	print(ratio)
	#f_loc = np.where(ratio > 1.03)
	print('earth cor')
	loc, leng = find_bad_reads(ratio)
	print('loc')
	recal(loc, file, leng)
	print('recall')
	os.remove(ima)
	os.remove(flt)
	print(os.getcwd())
	calwf3(file)
	time = gettime(file)
	header_update(file, time, loc)





def get_ratio(file):
	file_rhs = str(file + '[100:900,700:900]')
	file_lhs = str(file + '[100:900,50:250]')
	time,lhs = pstat.pstat(file_lhs,   units='rate', stat='midpt',diff=True,plot=False)
	time,rhs = pstat.pstat(file_rhs,  units='rate', stat='midpt',diff=True,plot=False)
	print(lhs)
	print(rhs)
	ratio=lhs/rhs
	return ratio
	

def find_bad_reads(ratio):
	f_loc = np.where(ratio > 1.05)
	print(f_loc)
	loc = f_loc[0]
	leng = len(ratio)
	return loc, leng


def recal(loc, file, leng):
	for r_loc in loc:
		if r_loc != (leng-1):
			read = r_loc + 1
			fits.setval(file,extver=read,extname='DQ',keyword='pixvalue',value=1024)


def gettime(file):
	flt = file[:-8]+'flt.fits'
	time = fits.getdata(flt, 5)
	median = np.median(time)
	return median


def header_update(file, median, loc):
	flt = file[:-8]+'flt.fits'
	hdu = fits.open(flt, mode='update')
	hdr = hdu[0].header
	exp = hdr['EXPTIME']
	hdr['OLDEXPT'] = exp
	hdr['EXPTIME'] = median
	location = loc + 1
	hdr.set('bad_read', str(location))
	hdu.close()


#def main():
#	file='/grp/hst/wfc3v/hkurtz/sky_flats/sin_test/icqtbbbxq_raw.fits'
#	e_lim_cor(file)
#main()






















