#file:blob_stack.py
#Description: Stacks the MDI files based on blob apearances one full image 
#per appearance. This can only be used if the 512 flag in not used when creating the mdi files.
#Date: Oct. 25, 2018
#author: Heather Kurtz

from astropy.io import ascii
from astropy.io import fits
import numpy as np
import numpy.ma as ma
import glob
import os

filters = ['F105W', 'F110W', 'F125W', 'F140W',  'F160W', 'F098M' ]

#'F132N', 'F139M', 'F126N', 'F153M', 'F127M','F128N', 'F164N'
#'F098M', 'F130N', 'F105W', 'F110W', 'F125W', 'F140W',  'F160W',

def Stack(data_array_1, signal_arr, time_arr, sum_arr, wht_arr, wavg_arr, mask, error):				#???
	im_median = np.nanmedian(data_array_1,axis=0)
	#image_median = im_median * mask
	im_mean   = np.nanmean(data_array_1,  axis=0)
	#image_mean   = im_mean * mask
	im_std    = np.nanstd(data_array_1,   axis=0)
	#image_std    = im_std * mask
	im_error    = np.sqrt(sum(error**2))/np.sqrt((len(error))) #np.nanmean(error,   axis=0)
	#image_error   = im_error * mask
	#im_min    = np.nanmin(data_array_1,   axis=0)
	#image_min    = im_min * mask
	im_signal = np.nansum(signal_arr, axis=0)
	#image_signal = im_signal * mask
	im_sum    = np.nansum(sum_arr,    axis=0)
	#image_sum    = im_sum * mask
	im_time   = np.nansum(time_arr,   axis=0)														
	#image_time   = im_time * mask														
	im_wht    = np.nansum(wht_arr,    axis=0)										#???											
	#image_wht    = im_wht * mask									#???											
	im_wavg   = np.nansum(wavg_arr,   axis=0)										#???											
	#image_wavg   = im_wavg * mask								#???											

	#return (image_mean, image_median, image_std, image_signal, image_time, image_sum, image_wht, image_wavg,image_error)#image_min
	return (im_mean, im_median, im_std, im_signal, im_time, im_sum, im_wht, im_wavg,im_error)#image_min
	#return (im_median)

def get_expstar(file,mjd_list):
	hdr = fits.getheader(file,0)
	mjd = hdr['EXPSTART']
	mjd_list.append(mjd) 
	return(mjd_list)


def read_in_bpix():
	filt_table = ascii.read('/user/holszewski/IR_flats/bpixtab_summary.txt', data_start=1, delimiter=',')
	bpixtab = filt_table['col1']
	usaftermjd = filt_table['col3']
	return (usaftermjd,bpixtab)


def main():

	#gets current directory
	current = os.getcwd()
	#base_path = '/user/hkurtz/IR_flats/test_f098' /grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_no512_run/
	base_path = '/grp/hst/wfc3v/hkurtz/sky_flats/Mar4_noflat_run/'
	#shange directory to the base path
	os.chdir(base_path)
	
	for fil in filters:
		print('working on ', fil)
		fil_path= os.path.join(base_path,fil)
		os.chdir(fil_path)
		#list of files
		list_of_files= glob.glob('*mdi.fits')
		hdr = fits.getheader(list_of_files[0], 1)
		nx = hdr['NAXIS1']
		ny = hdr['NAXIS2']
		#print(len(list_of_files))

#		for file in list_of_files:
#			hdr = fits.getheader(file,0)
#			mjd = hdr('EXPSTART')
#			if bpix_data[0] < mjd:
#				bpix_data0.append(file)
#			if bpix_data[1] < mjd:
#				bpix_data1.append(file)
#
#
		bpix_date_list, bpix = read_in_bpix()
		mjd_list = []
		date_dict = {}
		bpix_dict = {}
		for file in list_of_files:
			mjd_list = get_expstar(file,mjd_list)

		mjd_arr = np.array(mjd_list)

		files_arr = np.array(list_of_files)
		for i in range(len(bpix_date_list)):
			#this masks images with an observation date less than the date of the blob aprearances and stores them in a dictionary
			#print(bpix_date_list[i])
			date_mask = [bpix_date_list[i] < mjd_arr]
			date_dict[bpix_date_list[i]] = files_arr[date_mask]
			bpix_dict[bpix_date_list[i]] = bpix[i]

		for key in date_dict:
			# creats 3D array of size of the number of images in dictionary

			nf = len(date_dict[key])
			data_array = np.empty((nf, ny, nx), dtype=float)
			signal_arr = np.empty((nf, ny, nx), dtype=float)
			time_arr   = np.empty((nf, ny, nx), dtype=float)				
			sum_arr    = np.empty((nf, ny, nx), dtype=float)				
			err_arr    = np.empty((nf, ny, nx), dtype=float)									#???
			wht_arr    = np.empty((nf, ny, nx), dtype=float)									#???
			wavg_arr   = np.empty((nf, ny, nx), dtype=float)
			#print("array init")
			for i, f in enumerate(date_dict[key]):
				#reads in the data
				data=fits.getdata(f, 1)
				hdr = fits.getheader(f, 0)
				image = hdr['ROOTNAME']
				norm = hdr['NORM']
				exptime = hdr['EXPTIME']
				flag = hdr['EXPFLAG']
				#print('got header data')
				
       	
				if flag != 'INDETERMINATE':
					if exptime > 300:
						if "OLDEXPT" in hdr:
							old=hdr['OLDEXPT']
							diff=old-exptime
							#print('test data for outliers')
							if diff < 10.0:
								#print('***Good data:   ----- ',image, exptime,flag)
								data_array[i, :, :] = data
								signal_arr[i, :, :] = data * 0.0 + (norm * exptime)					##Corrected from Heather's v3
								time_arr[i, :, :]   = data * 0.0 + exptime				
								#sum_arr[i, :, :]    = data * 0.0 + 1		
								#err_arr[i, :, :]	= 1/(data * 0.0 + exptime*0.05 + 20.0**2 + norm*exptime)**0.5 	#???					
								#wht_arr[i, :, :]	= (signal_arr[i,:,:]*err_arr[i,:,:])								#???					
								#wavg_arr[i, :, :]	= data * (signal_arr[i,:,:]*err_arr[i,:,:])						#???
							else:																
								#print('***LIMB Rejected: --- ',image, old,exptime,diff)
								data_array[i, :, :] = np.nan
								signal_arr[i, :, :] = np.nan 
								time_arr[i, :, :]   = np.nan 				
								sum_arr[i, :, :]    = np.nan 
								err_arr[i, :, :]	= np.nan														 							
								wht_arr[i, :, :]	= np.nan																					
								wavg_arr[i, :, :]	= np.nan
						else:																
							#print('***good no header:   ----- ',image, exptime,flag)
							data_array[i, :, :] = data
							signal_arr[i, :, :] = data * 0.0 + (norm * exptime)					##Corrected from Heather's v3
							time_arr[i, :, :]   = data * 0.0 + exptime				
							#sum_arr[i, :, :]    = data * 0.0 + 1		
							#err_arr[i, :, :]	= 1/(data * 0.0 + exptime*0.05 + 20.0**2 + norm*exptime)**0.5 	#???					
							#wht_arr[i, :, :]	= (signal_arr[i,:,:]*err_arr[i,:,:])								#???					
							#wavg_arr[i, :, :]	= data * (signal_arr[i,:,:]*err_arr[i,:,:])	
					else:
						#print('***Bad EXPTIME: ----- ',image, exptime,flag)
						data_array[i, :, :] = np.nan
						signal_arr[i, :, :] = np.nan 
						time_arr[i, :, :]   = np.nan 				
						sum_arr[i, :, :]    = np.nan 
						err_arr[i, :, :]	= np.nan														 							
						wht_arr[i, :, :]	= np.nan																					
						wavg_arr[i, :, :]	= np.nan																			
				else:
					#print('***Bad EXPFLAG: ----- ',image, exptime,flag)
					data_array[i, :, :] = np.nan
					signal_arr[i, :, :] = np.nan 
					time_arr[i, :, :]   = np.nan 				
					sum_arr[i, :, :]    = np.nan 
					err_arr[i, :, :]	= np.nan														 							
					wht_arr[i, :, :]	= np.nan																					
					wavg_arr[i, :, :]	= np.nan		
							
					
				


			#run_path = os.getcwd()
			#mask_path = '/user/holszewski/IR_flats/'
			#os.chdir(mask_path)
			#print("data in array if nothing printed probelm")
			bpixtab = bpix_dict[key]
			blob_im_name = '/user/holszewski/IR_flats/' + bpixtab[:-11] + ".fits"
			blob_im_full = fits.getdata(blob_im_name)
			blob_im = blob_im_full[5:1019,5:1019]
			blob_im[blob_im < 511] = 0
			blob_im[blob_im > 511] = 1
			blob_im1full = fits.getdata('/user/holszewski/IR_flats/ground_blobs.fits')
			blob_im1 = blob_im1full[5:1019,5:1019]
			blob_im1[blob_im1 < 511] = 0
			blob_im1[blob_im1 > 511] = 1
			f_blob=blob_im1+blob_im
			f_blob[f_blob < 1] = 0
			f_blob[f_blob > 0] = 1
			#os.chdir(run_path)
			S_mean,S_median,S_std,S_signal,S_time,S_sum,S_wht,S_wavg,S_error=Stack(data_array, signal_arr, time_arr, sum_arr, wht_arr, wavg_arr, f_blob, err_arr)#S_min
			#S_median = Stack(data_array, signal_arr, time_arr, sum_arr, wht_arr, wavg_arr, f_blob, err_arr)
			#writing out the files
			n_mean = fil + '_' + str(key) + '_' + 'blob_mean_v2.fits'
			n_median = fil + '_' + str(key) + '_' + 'blob_median_v2.fits'
			n_sum = fil + '_' + str(key) + '_' + 'blob_sum_v2.fits'
			#n_min = fil + '_' + str(key) + '_' + 'blob_min_v2.fits'
			n_std = fil + '_' + str(key) + '_' + 'blob_std_v2.fits'
			n_signal = fil + '_' + str(key) + '_' + 'blob_signal_v2.fits'
			n_time = fil + '_' + str(key) + '_' + 'blob_time_v2.fits'
			n_error = fil + '_' + str(key) + '_' + 'blob_error_squred_v2.fits'
			fits.writeto(n_mean, S_mean,overwrite=True)
			fits.writeto(n_median, S_median,overwrite=True)
			fits.writeto(n_sum, S_sum,overwrite=True)
			fits.writeto(n_std, S_std,overwrite=True)
			#fits.writeto(n_min, S_min,overwrite=True)
			fits.writeto(n_signal, S_signal,overwrite=True)
			fits.writeto(n_time,   S_time,  overwrite=True)
			fits.writeto(n_error,   S_error,  overwrite=True)
		
main()































