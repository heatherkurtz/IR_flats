#Filename: stack_v3.py
#Description: This creates a stack of the masked, normalized frames 
#Date: August 12, 2018
#Author: Heather Kurtz

from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

filters = [ 'F098M', 'F105W', 'F110W', 'F125W', 'F140W', 'F160W' ]
#filters = [ 'F160W' ]

def Stack(data_array_1, signal_arr, time_arr, sum_arr, wht_arr, wavg_arr):				#???
	image_median = np.nanmedian(data_array_1,axis=0)
	image_mean   = np.nanmean(data_array_1,  axis=0)
	image_std    = np.nanstd(data_array_1,   axis=0)
	image_min    = np.nanmin(data_array_1,   axis=0)
	image_signal = np.nansum(signal_arr, axis=0)
	image_sum    = np.nansum(sum_arr,    axis=0)
	image_time   = np.nansum(time_arr,   axis=0)														
	image_wht    = np.nansum(wht_arr,    axis=0)										#???											
	image_wavg   = np.nansum(wavg_arr,   axis=0)										#???											

	return (image_mean, image_median, image_std, image_min, image_signal, image_time, image_sum, image_wht, image_wavg)			#???

def main():
	#gets current directory
	current = os.getcwd()
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/nodq/'
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/no512/'				 #Dec		(non-parallels)
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/el_he_105_140_run/'	 #March
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/apr_run/'				 #April
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/June10_run/'			 #June		bestrefs
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/test_new_files_run/'   #July  Noflat  (contains limb)
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/July3_noflat_asn_run/' #July  Flat  (contains limb), contains helium
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Aug6_run/'  		     #Aug Flat, limb corrected, helium corrected
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Aug6_noflat_run/'  		     #Aug Noflat, limb corrected, helium corrected
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Aug6_no512_run/'  		     #Aug Flat, limb corrected, helium corrected, blob flag off
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Aug6_no512_nolimb_run/'  		     #Aug Flat, contains limb, helium corrected, blob flag off

	base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Aug6_no512_noflat_run/'  		     #Aug Noflat, limb corrected, helium corrected, , blob flag off

	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/aug21_test_seg_blobs_run/' #Aug Flat, limb corrected, helium corrected testing blobs in seg
	#base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/aug26_test_flag512onlyinseg_run/'  #Aug Flat, limb corrected, helium corrected 512 flaged before seg
	#change directory to the base path
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
		time_arr   = np.empty((nf, ny, nx), dtype=float)				
		sum_arr    = np.empty((nf, ny, nx), dtype=float)				
		err_arr    = np.empty((nf, ny, nx), dtype=float)									#???
		wht_arr    = np.empty((nf, ny, nx), dtype=float)									#???
		wavg_arr   = np.empty((nf, ny, nx), dtype=float)									#???
		
		#read in the data and the DQs from the .fits for both extensions
		for i , f in enumerate(list_of_files):
		
			data_1=fits.getdata(f, 1)
			hdr = fits.getheader(f, 0)
			norm = hdr['NORM']
			exptime = hdr['EXPTIME']
			flag = hdr['EXPFLAG']
			if flag != 'INDETERMINATE':
				if exptime > 300:
					data_array[i, :, :] = data_1
					signal_arr[i, :, :] = data_1 * norm * exptime
					time_arr[i, :, :]   = data_1 * 0.0 + exptime				
					sum_arr[i, :, :]    = data_1 * 0.0 + 1		
					err_arr[i, :, :]	= 1/(data_1 * 0.0 + exptime*0.05 + 20.0**2 + norm*exptime)**0.5 		#???					
					wht_arr[i, :, :]	= (signal_arr[i,:,:]*err_arr[i,:,:])									#???					
					wavg_arr[i, :, :]	= data_1 * (signal_arr[i,:,:]*err_arr[i,:,:])
			else:
				print('short exposre or bad expflag')
				data_array[i, :, :] = np.nan
				signal_arr[i, :, :] = np.nan 
				time_arr[i, :, :]   = np.nan 				
				sum_arr[i, :, :]    = np.nan		
				err_arr[i, :, :]	= np.nan 		#???					
				wht_arr[i, :, :]	= np.nan									#???					
				wavg_arr[i, :, :]	= np.nan							#???			
					
			#test1=signal_arr[i,0,0]    #This is the lower left corner of the image
			#test2=   err_arr[i,0,0]
			#test3=   wht_arr[i,0,0]
			#print(norm,exptime,test1,test2,test3)											

#Norm				Exptime			Signal				Err					Wht=Sig/Err			
#														2x more error		4x lower weight
#0.8266774415969849 302.9384765625   250.432403564453 0.0387614555656909   9.707124482972745
#3.2508740425109860 652.9377441406  2122.618408203125 0.0197825375944376  41.990778458923614


# This part checks for exptime.  See /user/hkurtz/IR_flats/stack.py
			#data_1=fits.getdata(f, 1)
			#hdr = fits.getheader(f, 0)
			#norm = hdr['NORM']
			#exptime = hdr['EXPTIME']
			#if exptime > 300:
			#		data_array[i, :, :] = data_1
			#		signal_arr[i, :, :] = data_1 * 0.0 + (norm * exptime)
			#else:
			#		print('short exposre')
			#		data_array[i, :, :] = 0
			#		signal_arr[i, :, :] = 0 * norm * exptime

#This is a repeat that checks for expflag
            #data_1=fits.getdata(f, 1)
            #hdr = fits.getheader(f, 0)
            #norm = hdr['NORM']
            #exptime = hdr['EXPTIME']
            #flag = hdr['EXPFLAG']
            #if flag != 'indeterminate':
            #        if exptime > 300:
            #                data_array[i, :, :] = data_1
            #                signal_arr[i, :, :] = data_1 * norm * exptime		
            #        else:
            #                print('short exposre or bad expflag')
            #                data_array[i, :, :] = np.nan
            #                signal_arr[i, :, :] = np.nan * norm * exptime

			
		S_mean,S_median,S_std,S_min,S_signal,S_time,S_sum,S_wht,S_wavg=Stack(data_array, signal_arr, time_arr, sum_arr, wht_arr, wavg_arr)		

		n_mean   = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_mean.fits'
		n_median = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_median.fits'
		n_min    = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_min.fits'
		n_std    = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_std.fits'
		n_signal = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_signal.fits'
		n_sum    = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_noflat_sum.fits'
		n_time   = '/grp/hst/wfc3v/hkurtz/sky_flats/stacks/' + fil + '_aug6_no512_nofalt_time.fits'		
		#n_wht    = '/grp/hst/wfc3c/mack/ir_flat/PIRZKAL_2017/test_stack/' + fil + '_jul_wht.fits'			#???
		#n_wavg   = '/grp/hst/wfc3c/mack/ir_flat/PIRZKAL_2017/test_stack/' + fil + '_jul_wavg.fits'			#???

		fits.writeto(n_mean,   S_mean,  overwrite=True)
		fits.writeto(n_median, S_median,overwrite=True)
		fits.writeto(n_std,    S_std,   overwrite=True)
		fits.writeto(n_min,    S_min,   overwrite=True)
		fits.writeto(n_signal, S_signal,overwrite=True)
		fits.writeto(n_sum,    S_sum,   overwrite=True)
		fits.writeto(n_time,   S_time,  overwrite=True)					
		#fits.writeto(n_wht,    S_wht,   overwrite=True)														#???
		#fits.writeto(n_wavg,   S_wavg,  overwrite=True)														#???
	
main()
