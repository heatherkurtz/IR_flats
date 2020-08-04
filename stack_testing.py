#Filename: stack.py
#Description: This creates a mask of all the cosmic rays in the DQ and applies it to the image 
#Date: August 12, 2016
#Author: Heather Kurtz

from astropy.io import fits
from PyAstronomy import pyasl
import matplotlib.pylab as plt
import numpy as np
import numpy.ma as ma
import glob
import os

#filters = [ 'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W' ]
filters = ['early_140','ben_140']
#'F132N', 'F139M', 'F126N', 'F153M', 'F127M','F128N', 'F164N'
#'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W',

def Stack(data_array_1, signal_arr):
	image_median = np.nanmedian(data_array_1, axis=0)
	image_mean = np.nanmean(data_array_1, axis=0)
	image_sum = np.nansum(data_array_1, axis=0)
	image_signal = np.nansum(signal_arr, axis=0)
	image_std = np.nanstd(data_array_1, axis=0)
	image_min = np.nanmin(data_array_1, axis=0)

	return (image_mean, image_median, image_sum, image_std, image_min, image_signal)


#def number_stack(data_array_1):
	#get number of good values per pixel


def pixels(ar,per,nx,ny,start,end):#

    mean_im = np.zeros((ny, nx), dtype=float)
    median_im = np.zeros((ny, nx), dtype=float)
    sum_im = np.zeros((ny, nx), dtype=float)
    std_im = np.zeros((ny, nx), dtype=float)
    min_im = np.zeros((ny, nx), dtype=float)

    for i in range(len(ar[0,start:end,start:end])):
        for j in range(len(ar[0,0,start:end])):
            #pix = (ar[:,i,j])
            ni = i+start
            print(ni)
            nj = j+start
            print(nj)
            pix = (ar[:,ni,nj])
            print(len(pix))
            ok = np.isfinite(pix)
            pix = pix[ok]
            print(len(pix))
            file = 'before_pix_' + str(ni) + '_' + str(nj) + 'test_1.txt'
            np.savetxt(file, pix)
            bad = np.isnan(pix).sum()
            len_all = len(pix)
            val = int((len_all - bad)/2)
            if val < 2:
                print('bad_pixel')
            else:
                mask = np.isfinite(pix)
                pix = pix[mask]
                r = pyasl.generalizedESD(pix, val, per, fullOutput=True)
                print(r[1])
            plt.plot(pix, 'b.')
            for m in range(r[0]):
                plt.plot(r[1][m], pix[r[1][m]], 'r.')
            fname = 'pix_plot_' + str(ni) + '_' + str(nj) + 'test_25.jpg'
            plt.savefig(fname, dpi=300)
            plt.clf()
            pix[r[1][:]] = np.nan
            fileA = 'after_pix_' + str(ni) + '_' + str(nj) + 'test_25.txt'
            np.savetxt(fileA, pix)
            try:
            	mean_im[ni,nj]=np.nanmean(pix)
            	median_im[ni,nj]=np.nanmedian(pix)
            	sum_im[ni,nj]=np.nansum(pix)
            	std_im[ni,nj]=np.nanstd(pix)
            	min_im[ni,nj]=np.nanmin(pix)
            except ValueError: 
            	print('there is a value error', ni,nj)
    return(mean_im,median_im,sum_im,std_im,min_im)
            



def main():

	#gets current directory
	current = os.getcwd()
	#base_path = '/user/hkurtz/IR_flats/test_f098'
	base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/for_Ben/'
	#shange directory to the base path
	os.chdir(base_path)
	
	for fil in filters:
		fil_path= os.path.join(base_path,fil)
		os.chdir(fil_path)
		#list of files
		list_of_files= glob.glob('*mdi.fits')
		
		#gets the file size
		hdr = fits.getheader(list_of_files[0], 1)
		nx = hdr['NAXIS1']
		ny = hdr['NAXIS2']
		nf = len(list_of_files)
		set_data=fits.getdata(list_of_files[0], 1)
		#makes empty array to be filled with data
		data_array = np.empty((nf, ny, nx), dtype=float)
		signal_arr = np.empty((nf, ny, nx), dtype=float)
		#read in the data and the DQs from the .fits for both extensions
		for i , f in enumerate(list_of_files):
		
			data_1=fits.getdata(f, 1)
			hdr = fits.getheader(f, 0)
			norm = hdr['NORM']
			exptime = hdr['EXPTIME']
			data_array[i, :, :] = data_1
			signal_arr[i, :, :] = data_1 * norm * exptime
		start = 250
		end = 270
		per = 0.25
		S_mean,S_median,S_sum,S_std,S_min=pixels(data_array,per,nx,ny,start,end)
		#S_mean,S_median,S_sum,S_std,S_min, S_signal=Stack(data_array, signal_arr)
		n_mean = fil + '_mean_esd_test5nd_25.fits'
		n_median = fil + '_median_esd_test5nd_25.fits'
		n_sum = fil + '_sum_esd_test5nd_25.fits'
		n_min = fil + '_min_esd_test5nd_25.fits'
		n_std = fil + '_std_esd_test5nd_25.fits'
		n_signal = fil + '_signal_esd_test5nd_25.fits'
		fits.writeto(n_mean, S_mean,overwrite=True)
		fits.writeto(n_median, S_median,overwrite=True)
		fits.writeto(n_sum, S_sum,overwrite=True)
		fits.writeto(n_std, S_std,overwrite=True)
		fits.writeto(n_min, S_min,overwrite=True)
		#fits.writeto(n_signal, S_signal,overwrite=True)

	
main()




