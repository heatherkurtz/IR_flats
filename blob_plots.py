#file: blob_plots.py
#start with plot of blob 25

import matplotlib.pyplot as plt
import csv
import pandas as pd 
import glob
from astropy.io import ascii 
import numpy as np 
import os

#098,105,110,125,140,160
#lam = [9864,10552,11534,12486,13923,15369]
dict_lam = {'F098M':9864,'F105W':10552,'F110W':11534,'F125W':12486,'F140W':13923,'F160W':15369}
plot_path = '/grp/hst/wfc3v/hkurtz/sky_flats/mar4_noflat_plots/'

#def read_in_bpix():
#	filt_table = ascii.read('/user/holszewski/IR_flats/blob_summary_color.txt', data_start=0, names= ['n','x','y','radius','color','flux','mjd','flag'], delimiter=',')
#	#print(filt_table)
#	blobn = filt_table['n']
#	blobx = filt_table['x']
#	bloby = filt_table['y']
#	blobrad = filt_table['radius']
#	usaftermjd = filt_table['mjd']
#	flag =  filt_table['flag']
#	flux = filt_table['flux']
#	return (usaftermjd,blobx,bloby,blobrad,blobn,flag,flux)
#

def read_blob_sum():
	filt_table = ascii.read('/user/holszewski/IR_flats/blob_sumry_wide4.csv', data_start=1, 
							names= ['num','x','y','radius','flux','mjd',"app_win",'flag','mean140','mean160','mean110','mean125','mean105','mean098',
							'median140','median160','median110','median125','median105','median098'], delimiter=',')
	#print(filt_table)
	blobn = filt_table['num']
	blobx = filt_table['x']
	bloby = filt_table['y']
	blobrad = filt_table['radius']
	usaftermjd = filt_table['mjd']
	flag =  filt_table['flag']
	flux = filt_table['flux']
	f140 = filt_table['mean140']
	f160 = filt_table['mean160']
	f110 = filt_table['mean110']
	f125 = filt_table['mean125']
	f105 = filt_table['mean105']
	f098 = filt_table['mean098']
	md140 = filt_table['median140']
	md160 = filt_table['median160']
	md110 = filt_table['median110']
	md125 = filt_table['median125']
	md105 = filt_table['median105']
	md098 = filt_table['median098']
	return (usaftermjd,blobx,bloby,blobrad,blobn,flag,flux,f140,f160,f110,f125,f105,f098,md140,md160,md110,md125,md105,md098)


#def read_in_bpix098():
#	filt_table = ascii.read('/user/holszewski/IR_flats/blob_summary_new_copy.txt', data_start=2, names= ['n','x','y','radius','flux','mjd','window'],delimiter=' ' )
#	#print(filt_table)
#	blobn = filt_table['n']
#	blobx = filt_table['x']
#	bloby = filt_table['y']
#	blobrad = filt_table['radius']
#	usaftermjd = filt_table['mjd']
#	#flag =  filt_table['Flag']
#	return (usaftermjd,blobx,bloby,blobrad,blobn, flag)
#

def main():
	print('start')
	#mjd,blobx,bloby,blobrad,blobn,flags,flux = read_in_bpix()
	#mjd,blobx,bloby,blobrad,blobn,flags = read_in_bpix098()
	mjd,blobx,bloby,blobrad,blobn,flags,flux,f140,f160,f110,f125,f105,f098,md140,md160,md110,md125,md105,md098 = read_blob_sum()

	print('read bpix_summery')
	mean_list = glob.glob('F*blob_color_mean_final_v3_noT0.csv')
	sum_list = glob.glob('F*blob_color_sum_final.csv')
	min_list = glob.glob('F*blob_color_min_final_v3_noT0.csv')
	plt.clf()

	for i in range(len(mean_list)):
		mf = mean_list[i]
		filt = mf[:5]
		lam = dict_lam[filt]
		df_mean = pd.read_csv(mean_list[i],names=['id', 'mean'])
		df_sum = pd.read_csv(sum_list[i],names=['id', 'sum'])
		mean_array = np.array(df_mean['mean'])
		ar = np.ones_like(mean_array)
		lam_ar = ar*lam
		if filt == 'F098M':
			f98means = mean_array/f098
			plt.plot(lam_ar, f98means, '.',  label=mf[:5])
		if filt == 'F105W':
			fm105eans = mean_array/f105
			plt.plot(lam_ar, fm105eans, '.',  label=mf[:5])
		if filt == 'F125W':
			f125means = mean_array/f125
			plt.plot(lam_ar, f125means, '.',  label=mf[:5])
		if filt == 'F110W':
			f110means = mean_array/f110
			plt.plot(lam_ar, f110means, '.',  label=mf[:5])
		if filt == 'F140W':
			f140means = mean_array/f140
			plt.plot(lam_ar, f140means, '.',  label=mf[:5])
		if filt == 'F160W':
			f160means = mean_array/f160
			plt.plot(lam_ar, f160means, '.',  label=mf[:5])

		#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
		#plt.plot(lam_ar, df_mean['mean'], '.', )# label=mf[:5])
	plt.xlabel('Pivit Wavelength')
	#plt.xlim(0,0.2)	
	plt.ylabel('Mean')
	plt.title(' Blob Mean per Wavelength')	
	#plt.legend(prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Mean_per_wavelength_noT0_v4_no_norm.png'))
	plt.clf()


	for i in range(len(mean_list)):
		mf = mean_list[i]
		filt = mf[:5]
		lam = dict_lam[filt]
		df_mean = pd.read_csv(mean_list[i],names=['id', 'mean'])
		df_sum = pd.read_csv(sum_list[i],names=['id', 'sum'])
		df_min = pd.read_csv(min_list[i],names=['id', 'min'])
		min_array = np.array(df_min['min'])
		ar = np.ones_like(min_array)
		lam_ar = ar*lam
		if filt == 'F098M':
			f98min = min_array/f098
			plt.plot(lam_ar, f98min, '.',  label=mf[:5])
		if filt == 'F105W':
			fm105in = min_array/f105
			plt.plot(lam_ar, fm105in, '.',  label=mf[:5])
		if filt == 'F125W':
			f125min = min_array/f125
			plt.plot(lam_ar, f125min, '.',  label=mf[:5])
		if filt == 'F110W':
			f110min = min_array/f110
			plt.plot(lam_ar, f110min, '.',  label=mf[:5])
		if filt == 'F140W':
			f140min = min_array/f140
			plt.plot(lam_ar, f140min, '.',  label=mf[:5])
		if filt == 'F160W':
			f160min = min_array/f160
			plt.plot(lam_ar, f160min, '.',  label=mf[:5])

		#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
		#plt.plot(lam_ar, df_min['min'], '.',)# label=mf[:5])
	plt.xlabel('Pivit Wavelength')
	#plt.xlim(0,0.2)	
	plt.ylabel('Minimum')
	plt.title(' Blob Minimum per Wavelength')	
	#plt.legend(prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_noT0_v4_no_norm.png'))
	plt.clf()



	#filter ratios
	f098m = 'F098M_blob_color_mean_final_noT0.csv'
	f105w = 'F105W_blob_color_mean_final_noT0.csv'
	f110w = 'F110W_blob_color_mean_final_noT0.csv'
	f125w = 'F125W_blob_color_mean_final_noT0.csv'
	f140w = 'F140W_blob_color_mean_final_noT0.csv'
	f160w = 'F160W_blob_color_mean_final_noT0.csv'


	df_098 = pd.read_csv(f098m,names=['id', 'mean'])
	df_105 = pd.read_csv(f105w,names=['id', 'mean'])
	df_110 = pd.read_csv(f110w,names=['id', 'mean'])
	df_125 = pd.read_csv(f125w,names=['id', 'mean'])
	df_140 = pd.read_csv(f140w,names=['id', 'mean'])
	df_160 = pd.read_csv(f160w,names=['id', 'mean'])

	r125_160=(df_125['mean']/f125)/(df_160['mean']/f160)
	r105_125=(df_105['mean']/f105)/(df_125['mean']/f125)
	r098_125=(df_098['mean']/f098)/(df_125['mean']/f125)
	r110_125=(df_110['mean']/f110)/(df_125['mean']/f125)
	#cm = plt.cm.get_cmap('CMRmap')
	cmap=plt.cm.get_cmap('rainbow', 20)
	plt.rcParams['axes.facecolor'] = 'lightgray'

	z=flux#/(max(flux))#df_160['mean']


	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(r125_160,r105_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('105/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_105_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(r125_160,r098_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_098_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	sc = plt.scatter(r125_160,r110_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('110/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_110_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	
	mean125 = df_125['mean']/f125
	mean160 = df_160['mean']/f160
	mean098 = df_098['mean']/f098
	mean105 = df_105['mean']/f105
	mean110 = df_110['mean']/f110
	mean140 = df_140['mean']/f140
	test = df_125['mean']
	print('norm =', mean125[0], 'ori =', test[0])


	yes125 = []
	no125 = []
	yes160 = []
	no160 = []
	yes098 = []
	no098 = []
	yes105 = []
	no105 = []
	yes110 = []
	no110 = []
	yes140 = []
	no140 = []
	for i in range(len(mean125)):
		#print(flags[i])
		if flags[i] == 'Yes':
			yes125.append(mean125[i])
			yes160.append(mean160[i])
			yes098.append(mean098[i])
			yes105.append(mean105[i])
			yes110.append(mean110[i])
			yes140.append(mean140[i])
		else:
			no125.append(mean125[i])
			no160.append(mean160[i])
			no098.append(mean098[i])
			no105.append(mean105[i])
			no110.append(mean110[i])
			no140.append(mean140[i])
	aryes125 =np.array(yes125)
	arno125 =np.array(no125)
	aryes160=np.array(yes160)
	arno160 =np.array(no160)
	aryes098=np.array(yes098)
	arno098 =np.array(no098)
	aryes105 =np.array(yes105)
	arno105 =np.array(no105)
	aryes110=np.array(yes110)
	arno110 =np.array(no110)
	aryes140=np.array(yes140)
	arno140 =np.array(no140)
	#change plot opasity to see overlapping point
	#add plot with blob stength/flux(mean) as thierd axis

	yesr125_160=aryes125/aryes160
	#print(len(yesr125_160))
	nor125_160=arno125/arno160
	#print(len(nor125_160))
	yesr098_125=aryes098/aryes125
	nor098_125=arno098/arno125
	plt.rcParams['axes.facecolor'] = 'white'

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(yesr125_160,yesr098_125, '.', label='Flagged')
	plt.plot(nor098_125, nor125_160, 's', label='Not Flagged')
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Mean Colors With Flags')	
	plt.legend(prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_colors_flags_noT0_v4_no_norm.png'))
	plt.clf()


	#27,140,33,25,9,30,148
	blob25 = [mean098[26],  mean105[26], mean110[26], mean125[26], mean140[26], mean160[26]]
	blob100 = [mean098[1],mean105[1],mean110[1],mean125[1],mean140[1],mean160[1]]
	blob65 = [mean098[24],  mean105[24], mean110[24], mean125[24], mean140[24], mean160[24]]
	blob98 = [mean098[32],  mean105[32], mean110[32], mean125[32], mean140[32], mean160[32]]
	blob15 = [mean098[8],   mean105[8],  mean110[8],  mean125[8],  mean140[8],  mean160[8]]
	blob30 = [mean098[29],  mean105[29], mean110[29], mean125[29], mean140[29], mean160[29]]
	blob130 = [mean098[14],mean105[14],mean110[14],mean125[14],mean140[14],mean160[14]]

	lam_list = [9864,10552,11534,12486,13923,15369]




	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(lam_list, blob25, '.', label='blob27')
	plt.plot(lam_list, blob100, '.', label='blob2')
	plt.plot(lam_list, blob65, '.', label='blob25')
	plt.plot(lam_list, blob98, '.', label='blob33')
	plt.plot(lam_list, blob15, '.', label='blob9')
	plt.plot(lam_list, blob30, '.', label='blob30')
	plt.plot(lam_list, blob130, '.', label='blob15')
	plt.xlabel('Pivit Wavelength')
	plt.ylim(0.94,1.02)	
	plt.ylabel('Mean')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Bright_Blob_Mean_sub_set_noT0_v4_no_norm.png'))
	plt.clf()

	 

	print(len(aryes105))
	plt.plot(mean105, mean098, 'o',label='F098', markersize=4)
	plt.plot(mean105, mean160, '*',label='F160', markersize=4)
	plt.plot(mean105, mean110, 's',label='F110', markersize=4)
	plt.plot(mean105, mean125, 'd',label='F125', markersize=4)
	plt.plot(mean105, mean140, 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_all_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mean105, mean098 ,'o',label='F098', markersize=4)
	plt.plot(mean105, mean160 ,'*',label='F160', markersize=4)
	plt.plot(mean105, mean110 ,'s',label='F110', markersize=4)
	plt.plot(mean105, mean125 ,'d',label='F125', markersize=4)
	plt.plot(mean105, mean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_all_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot( mean105, mean098 ,'o',label='F098', markersize=4)
	plt.plot( mean105, mean160 ,'*',label='F160', markersize=4)
	plt.plot( mean105, mean110 ,'s',label='F110', markersize=4)
	plt.plot( mean105, mean125 ,'d',label='F125', markersize=4)
	plt.plot( mean105, mean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.9, 1.1 )
	plt.ylim(0.9, 1.1 )	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Small_Zoomed_Blob_Mean_vs_105_all_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	print(len(aryes105))
	plt.plot( aryes105, aryes098 , 'o',label='F098', markersize=4)
	plt.plot( aryes105, aryes160 , '*',label='F160', markersize=4)
	plt.plot( aryes105, aryes110 , 's',label='F110', markersize=4)
	plt.plot( aryes105, aryes125 , 'd',label='F125', markersize=4)
	plt.plot( aryes105, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_all_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(aryes105, aryes098, 'o',label='F098', markersize=4)
	plt.plot(aryes105, aryes160, '*',label='F160', markersize=4)
	plt.plot(aryes105, aryes110, 's',label='F110', markersize=4)
	plt.plot(aryes105, aryes125, 'd',label='F125', markersize=4)
	plt.plot(aryes105, aryes140, 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_all_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot(aryes105, aryes160 ,  '*',label='F160', markersize=4)
	plt.plot(aryes105, aryes110 ,  's',label='F110', markersize=4)
	plt.plot(aryes105, aryes125 ,  'd',label='F125', markersize=4)
	plt.plot(aryes105, aryes140 ,  'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot( aryes110, aryes160 , '*',label='F160', markersize=4)
	plt.plot( aryes110, aryes105 , 's',label='F105', markersize=4)
	plt.plot( aryes110, aryes125 , 'd',label='F125', markersize=4)
	plt.plot( aryes110, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_110_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot( aryes125, aryes160 ,'*',label='F160', markersize=4)
	plt.plot( aryes125, aryes110 ,'s',label='F110', markersize=4)
	plt.plot( aryes125, aryes105 ,'d',label='F105', markersize=4)
	plt.plot( aryes125, aryes140 ,'p',label='F140', markersize=4)
	plt.plot( aryes125, aryes125 ,'-',label='F125', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_125_flaged_noT0_v4_no_norm_line_final_v1.png'))
	plt.clf()


	plt.plot( aryes140, aryes160 , '*',label='F160', markersize=4)
	plt.plot( aryes140, aryes110 , 's',label='F110', markersize=4)
	plt.plot( aryes140, aryes125 , 'd',label='F125', markersize=4)
	plt.plot( aryes140, aryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_140_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	m105, b105 = np.polyfit(aryes160, aryes105, 1)
	slope105 = m105*aryes160+b105
	m110, b110 = np.polyfit(aryes160, aryes110, 1)
	slope110 = m110*aryes160+b110
	m140, b140 = np.polyfit(aryes160, aryes140, 1)
	slope140 = m140*aryes160+b140
	m125, b125 = np.polyfit(aryes160, aryes125, 1)
	slope125 = m125*aryes160+b125


	plt.plot(aryes160, aryes105 ,  'b*',label='F105', markersize=4)
	plt.plot(aryes160, slope105, 'b-',label='Fit F105')
	plt.plot(aryes160, aryes110 ,  'gs',label='F110', markersize=4)
	plt.plot(aryes160, slope110, 'g--',label='Fit F110')
	plt.plot(aryes160, aryes125 ,  'cd',label='F125', markersize=4)
	plt.plot(aryes160, slope125, 'c-.',label='Fit F125')
	plt.plot(aryes160, aryes140 ,  'yp',label='F140', markersize=4)
	plt.plot(aryes160, slope140, 'y:',label='Fit F140')
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_160_flaged_noT0_test_fit_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot( aryes105, aryes160 ,'*',label='F160', markersize=4)
	plt.plot( aryes105, aryes110 ,'s',label='F110', markersize=4)
	plt.plot( aryes105, aryes125 ,'d',label='F125', markersize=4)
	plt.plot( aryes105, aryes140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot( aryes110, aryes160 ,'*',label='F160', markersize=4)
	plt.plot( aryes110, aryes105 ,'s',label='F105', markersize=4)
	plt.plot( aryes110, aryes125 ,'d',label='F125', markersize=4)
	plt.plot( aryes110, aryes140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_110_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(aryes125, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes125, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes125, aryes105 , 'd',label='F105', markersize=4)
	plt.plot(aryes125, aryes125 , '-',label='F125', markersize=4)
	plt.plot(aryes125, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_125_flaged_noT0_v4_no_norm_line_test.png'))
	plt.clf()


	plt.plot( aryes140, aryes160 , '*',label='F160', markersize=4)
	plt.plot( aryes140, aryes110 , 's',label='F110', markersize=4)
	plt.plot( aryes140, aryes125 , 'd',label='F125', markersize=4)
	plt.plot( aryes140, aryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_140_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	m105, b105 = np.polyfit(aryes160, aryes105, 1)
	slope105 = m105*aryes160+b105
	m110, b110 = np.polyfit(aryes160, aryes110, 1)
	slope110 = m110*aryes160+b110
	m140, b140 = np.polyfit(aryes160, aryes140, 1)
	slope140 = m140*aryes160+b140
	m125, b125 = np.polyfit(aryes160, aryes125, 1)
	slope125 = m125*aryes160+b125


	plt.plot(aryes160, aryes105 ,  'b*',label='F105', markersize=4)
	plt.plot(aryes160, slope105, 'b-',label='Fit F105')
	plt.plot(aryes160, aryes110 ,  'gs',label='F110', markersize=4)
	plt.plot(aryes160, slope110, 'g--',label='Fit F110')
	plt.plot(aryes160, aryes125 ,  'cd',label='F125', markersize=4)
	plt.plot(aryes160, slope125, 'c-.',label='Fit F125')
	plt.plot(aryes160, aryes140 ,  'yp',label='F140', markersize=4)
	plt.plot(aryes160, slope140, 'y:',label='Fit F140')
	plt.xlabel('Mean of 160 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_160_flaged_noT0_fit_v4_no_norm.png'))
	plt.clf()


	plt.plot( mean110, mean160 , '*',label='F160', markersize=4)
	plt.plot( mean110, mean105 , 's',label='F105', markersize=4)
	plt.plot( mean110, mean125 , 'd',label='F125', markersize=4)
	plt.plot( mean110, mean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_110_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot( mean125, mean160 ,'*',label='F160', markersize=4)
	plt.plot( mean125, mean110 ,'s',label='F110', markersize=4)
	plt.plot( mean125, mean105 ,'d',label='F105', markersize=4)
	plt.plot( mean125, mean125 ,'-',label='F125', markersize=4)
	plt.plot( mean125, mean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_125_noT0_v4_no_norm_line_test.png'))
	plt.clf()


	plt.plot( mean140, mean160 , '*',label='F160', markersize=4)
	plt.plot( mean140, mean110 , 's',label='F110', markersize=4)
	plt.plot( mean140, mean125 , 'd',label='F125', markersize=4)
	plt.plot( mean140, mean105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_140_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mean160, mean105 ,  '*',label='F105', markersize=4)
	plt.plot(mean160, mean110 ,  's',label='F110', markersize=4)
	plt.plot(mean160, mean125 ,  'd',label='F125', markersize=4)
	plt.plot(mean160, mean140 ,  'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_160_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot( mean105, mean160 ,'*',label='F160', markersize=4)
	plt.plot( mean105, mean110 ,'s',label='F110', markersize=4)
	plt.plot( mean105, mean125 ,'d',label='F125', markersize=4)
	plt.plot( mean105, mean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot( mean110, mean160 ,'*',label='F160', markersize=4)
	plt.plot( mean110, mean105 ,'s',label='F105', markersize=4)
	plt.plot( mean110, mean125 ,'d',label='F125', markersize=4)
	plt.plot( mean110, mean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_110_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mean125, mean160 , '*',label='F160', markersize=4)
	plt.plot(mean125, mean110 , 's',label='F110', markersize=4)
	plt.plot(mean125, mean105 , 'd',label='F105', markersize=4)
	plt.plot(mean125, mean125 , '-',label='F125', markersize=4)
	plt.plot(mean125, mean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_125_noT0_v4_no_norm_line_test.png'))
	plt.clf()


	plt.plot( mean140, mean160 , '*',label='F160', markersize=4)
	plt.plot( mean140, mean110 , 's',label='F110', markersize=4)
	plt.plot( mean140, mean125 , 'd',label='F125', markersize=4)
	plt.plot( mean140, mean105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_140_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mean160, mean105 , '*',label='F105', markersize=4)
	plt.plot(mean160, mean110 , 's',label='F110', markersize=4)
	plt.plot(mean160, mean125 , 'd',label='F125', markersize=4)
	plt.plot(mean160, mean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_160_noT0_v4_no_norm.png'))
	plt.clf()


	#mins
	#filter ratiosblob_color_mean_final_noT0
	fmin098m = 'F098M_blob_color_min_final_noT0.csv'
	fmin105w = 'F105W_blob_color_min_final_noT0.csv'
	fmin110w = 'F110W_blob_color_min_final_noT0.csv'
	fmin125w = 'F125W_blob_color_min_final_noT0.csv'
	fmin140w = 'F140W_blob_color_min_final_noT0.csv'
	fmin160w = 'F160W_blob_color_min_final_noT0.csv'


	min_098 = pd.read_csv(fmin098m,names=['id', 'min'])
	min_105 = pd.read_csv(fmin105w,names=['id', 'min'])
	min_110 = pd.read_csv(fmin110w,names=['id', 'min'])
	min_125 = pd.read_csv(fmin125w,names=['id', 'min'])
	min_140 = pd.read_csv(fmin140w,names=['id', 'min'])
	min_160 = pd.read_csv(fmin160w,names=['id', 'min'])

	minr125_160=(min_125['min']/f125)/(min_160['min']/f160)
	minr105_125=(min_105['min']/f105)/(min_125['min']/f125)
	minr098_125=(min_098['min']/f098)/(min_125['min']/f125)
	minr110_125=(min_110['min']/f110)/(min_125['min']/f125)



	#cm = plt.cm.get_cmap('CMRmap')
	cmap=plt.cm.get_cmap('rainbow', 20)
	plt.rcParams['axes.facecolor'] = 'lightgray'

	z=flux#/(max(flux))#df_160['mean']


	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(minr125_160,minr105_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('105/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_105_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(minr125_160,minr098_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_098_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	sc = plt.scatter(minr125_160,minr110_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('110/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_110_125v125_160_noT0_v4_no_norm.png'))
	plt.clf()

	
	min125 = min_125['min']/f125
	min160 = min_160['min']/f160
	min098 = min_098['min']/f098
	min105 = min_105['min']/f105
	min110 = min_110['min']/f110
	min140 = min_140['min']/f140
	test = min_125['min']
	print('norm =', min125[0], 'ori =', test[0])


	yes125 = []
	no125 = []
	yes160 = []
	no160 = []
	yes098 = []
	no098 = []
	yes105 = []
	no105 = []
	yes110 = []
	no110 = []
	yes140 = []
	no140 = []
	for i in range(len(min125)):
		#print(flags[i])
		if flags[i] == 'Yes':
			yes125.append(min125[i])
			yes160.append(min160[i])
			yes098.append(min098[i])
			yes105.append(min105[i])
			yes110.append(min110[i])
			yes140.append(min140[i])
		else:
			no125.append(min125[i])
			no160.append(min160[i])
			no098.append(min098[i])
			no105.append(min105[i])
			no110.append(min110[i])
			no140.append(min140[i])
	aryes125 =np.array(yes125)
	arno125 =np.array(no125)
	aryes160=np.array(yes160)
	arno160 =np.array(no160)
	aryes098=np.array(yes098)
	arno098 =np.array(no098)
	aryes105 =np.array(yes105)
	arno105 =np.array(no105)
	aryes110=np.array(yes110)
	arno110 =np.array(no110)
	aryes140=np.array(yes140)
	arno140 =np.array(no140)
	#change plot opasity to see overlapping point
	#add plot with blob stength/flux(mean) as thierd axis

	yesr125_160=aryes125/aryes160
	#print(len(yesr125_160))
	nor125_160=arno125/arno160
	#print(len(nor125_160))
	yesr098_125=aryes098/aryes125
	nor098_125=arno098/arno125
	plt.rcParams['axes.facecolor'] = 'white'

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(yesr125_160,yesr098_125, '.', label='Flagged')
	plt.plot(nor098_125, nor125_160, 's', label='Not Flagged')
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Min Colors With Flags')	
	plt.legend()
	plt.savefig(os.path.join(plot_path,'Min_Blob_colors_flags_noT0_v4_no_norm.png'))
	plt.clf()


	#27,140,33,25,9,30,148
	blob25 = [min098[26], min105[26],  min110[26], min125[26], min140[26], min160[26]]
	blob100 =[min098[1],min105[1], min110[1],min125[1],min140[1],min160[1]]
	blob65 = [min098[24], min105[24],  min110[24], min125[24], min140[24], min160[24]]
	blob98 = [min098[32], min105[32],  min110[32], min125[32], min140[32], min160[32]]
	blob15 = [min098[8],  min105[8],   min110[8],  min125[8],  min140[8],  min160[8]]
	blob30 = [min098[29], min105[29],  min110[29], min125[29], min140[29], min160[29]]
	blob130 =[min098[14],min105[14], min110[14],min125[14],min140[14],min160[14]]

	lam_list = [9864,10552,11534,12486,13923,15369]




	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(lam_list, blob25, '.', label='blob27')
	plt.plot(lam_list, blob100, '.', label='blob2')
	plt.plot(lam_list, blob65, '.', label='blob25')
	plt.plot(lam_list, blob98, '.', label='blob33')
	plt.plot(lam_list, blob15, '.', label='blob9')
	plt.plot(lam_list, blob30, '.', label='blob30')
	plt.plot(lam_list, blob130, '.', label='blob15')
	plt.xlabel('Pivit Wavelength')
	#plt.xlim(0,0.2)	
	plt.ylabel('Min')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Bright_Blob_Min_sub_set_noT0_v4_no_norm.png'))
	plt.clf()

	 

	print(len(aryes105))
	plt.plot( aryes105, aryes098 , 'o',label='F098', markersize=4)
	plt.plot( aryes105, aryes160 , '*',label='F160', markersize=4)
	plt.plot( aryes105, aryes110 , 's',label='F110', markersize=4)
	plt.plot( aryes105, aryes125 , 'd',label='F125', markersize=4)
	plt.plot( aryes105, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Min in each filter')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_105_all_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot( aryes105, aryes098, 'o',label='F098', markersize=4)
	plt.plot( aryes105, aryes160, '*',label='F160', markersize=4)
	plt.plot( aryes105, aryes110, 's',label='F110', markersize=4)
	plt.plot( aryes105, aryes125, 'd',label='F125', markersize=4)
	plt.plot( aryes105, aryes140, 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )	
	plt.ylabel('Min in each filter')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best', prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_105_all_noT0_v4_no_norm.png'))
	plt.clf()



	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	print(len(aryes160),len(aryes105),len(aryes110),len(aryes125),len(aryes140))
	plt.plot( aryes105, aryes160 ,'*',label='F160', markersize=4)
	plt.plot( aryes105, aryes110 ,'s',label='F110', markersize=4)
	plt.plot( aryes105, aryes125 ,'d',label='F125', markersize=4)
	plt.plot( aryes105, aryes140 ,'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 105')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_105_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(aryes110, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes110, aryes105 , 's',label='F105', markersize=4)
	plt.plot(aryes110, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes110, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 110 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 110')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_110_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(aryes125, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes125, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes125, aryes105 , 'd',label='F105', markersize=4)
	plt.plot(aryes125, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 125 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 125')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_125_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(aryes140, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes140, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes140, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes140, aryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Min of 140 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 140')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_140_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(aryes160, aryes105 , '*',label='F105', markersize=4)
	plt.plot(aryes160, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes160, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes160, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 160')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_160_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot(aryes105, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes105, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes105, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes105, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 105')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_105_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(aryes110, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes110, aryes105 , 's',label='F105', markersize=4)
	plt.plot(aryes110, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes110, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 110 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 110')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_110_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(aryes125, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes125, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes125, aryes105 , 'd',label='F105', markersize=4)
	plt.plot(aryes125, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 125 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 125')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_125_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(aryes140, aryes160 , '*',label='F160', markersize=4)
	plt.plot(aryes140, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes140, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes140, aryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Min of 140 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 140')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_140_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(aryes160, aryes105 , '*',label='F105', markersize=4)
	plt.plot(aryes160, aryes110 , 's',label='F110', markersize=4)
	plt.plot(aryes160, aryes125 , 'd',label='F125', markersize=4)
	plt.plot(aryes160, aryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 160 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 160')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_160_noT0_v4_no_norm.png'))
	plt.clf()

	for i in range(len(mean_list)):
		mf = mean_list[i]
		filt = mf[:5]
		lam = dict_lam[filt]
		df_mean = pd.read_csv(mean_list[i],names=['id', 'mean'])
		df_sum = pd.read_csv(sum_list[i],names=['id', 'sum'])
		mean_array = np.array(df_mean['mean'])
		ar = np.ones_like(mean_array)
		lam_ar = ar*lam
		if filt == 'F098M':
			md98means = mean_array/md098
			plt.plot(lam_ar, md98means, '.',  label=mf[:5])
		if filt == 'F105W':
			mdm105eans = mean_array/md105
			plt.plot(lam_ar, mdm105eans, '.',  label=mf[:5])
		if filt == 'F125W':
			md125means = mean_array/md125
			plt.plot(lam_ar, md125means, '.',  label=mf[:5])
		if filt == 'F110W':
			md110means = mean_array/md110
			plt.plot(lam_ar, md110means, '.',  label=mf[:5])
		if filt == 'F140W':
			md140means = mean_array/md140
			plt.plot(lam_ar, md140means, '.',  label=mf[:5])
		if filt == 'F160W':
			md160means = mean_array/md160
			plt.plot(lam_ar, md160means, '.',  label=mf[:5])

		#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
		#plt.plot(lam_ar, df_mean['mean'], '.', )# label=mf[:5])
	plt.xlabel('Pivit Wavelength')
	plt.ylabel('Mean')
	plt.title(' Blob Mean per Wavelength')	
	plt.savefig(os.path.join(plot_path,'Blob_Mean_per_wavelength_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	for i in range(len(mean_list)):
		mf = mean_list[i]
		filt = mf[:5]
		lam = dict_lam[filt]
		df_mean = pd.read_csv(mean_list[i],names=['id', 'mean'])
		df_sum = pd.read_csv(sum_list[i],names=['id', 'sum'])
		df_min = pd.read_csv(min_list[i],names=['id', 'min'])
		min_array = np.array(df_min['min'])
		ar = np.ones_like(min_array)
		lam_ar = ar*lam
		if filt == 'F098M':
			md98min = min_array/md098
			plt.plot(lam_ar, md98min, '.',  label=mf[:5])
		if filt == 'F105W':
			mdm105in = min_array/md105
			plt.plot(lam_ar, mdm105in, '.',  label=mf[:5])
		if filt == 'F125W':
			md125min = min_array/md125
			plt.plot(lam_ar, md125min, '.',  label=mf[:5])
		if filt == 'F110W':
			md110min = min_array/md110
			plt.plot(lam_ar, md110min, '.',  label=mf[:5])
		if filt == 'F140W':
			md140min = min_array/md140
			plt.plot(lam_ar, md140min, '.',  label=mf[:5])
		if filt == 'F160W':
			md160min = min_array/md160
			plt.plot(lam_ar, md160min, '.',  label=mf[:5])

		#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
		#plt.plot(lam_ar, df_min['min'], '.',)# label=mf[:5])
	plt.xlabel('Pivit Wavelength')
	#plt.xlim(0,0.2)	
	plt.ylabel('Minimum')
	plt.title(' Blob Minimum per Wavelength')	
	#plt.legend(prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_md_norm_noT0_v4_no_norm.png'))
	plt.clf()



	#filter ratios
	f098m = 'F098M_blob_color_mean_final_noT0.csv'
	f105w = 'F105W_blob_color_mean_final_noT0.csv'
	f110w = 'F110W_blob_color_mean_final_noT0.csv'
	f125w = 'F125W_blob_color_mean_final_noT0.csv'
	f140w = 'F140W_blob_color_mean_final_noT0.csv'
	f160w = 'F160W_blob_color_mean_final_noT0.csv'


	df_098 = pd.read_csv(f098m,names=['id', 'mean'])
	df_105 = pd.read_csv(f105w,names=['id', 'mean'])
	df_110 = pd.read_csv(f110w,names=['id', 'mean'])
	df_125 = pd.read_csv(f125w,names=['id', 'mean'])
	df_140 = pd.read_csv(f140w,names=['id', 'mean'])
	df_160 = pd.read_csv(f160w,names=['id', 'mean'])

	md125_160=(df_125['mean']/md125)/(df_160['mean']/md160)
	md105_125=(df_105['mean']/md105)/(df_125['mean']/md125)
	md098_125=(df_098['mean']/md098)/(df_125['mean']/md125)
	md110_125=(df_110['mean']/md110)/(df_125['mean']/md125)
	#cm = plt.cm.get_cmap('CMRmap')
	cmap=plt.cm.get_cmap('rainbow', 20)
	plt.rcParams['axes.facecolor'] = 'lightgray'

	z=flux#/(max(flux))#df_160['mean']


	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(md125_160,md105_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('105/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_105_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(md125_160,md098_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_098_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	sc = plt.scatter(md125_160,md110_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('110/125')
	plt.title(' Blob Mean Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Mean_colors_110_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	
	mdmean125 = df_125['mean']/md125
	mdmean160 = df_160['mean']/md160
	mdmean098 = df_098['mean']/md098
	mdmean105 = df_105['mean']/md105
	mdmean110 = df_110['mean']/md110
	mdmean140 = df_140['mean']/md140
	test = df_125['mean']
	print('norm =', mdmean125[0], 'ori =', test[0])


	mdyes125 = []
	mdno125 = []
	mdyes160 = []
	mdno160 = []
	mdyes098 = []
	mdno098 = []
	mdyes105 = []
	mdno105 = []
	mdyes110 = []
	mdno110 = []
	mdyes140 = []
	mdno140 = []
	for i in range(len(mdmean125)):
		#print(flags[i])
		if flags[i] == 'Yes':
			mdyes125.append(mdmean125[i])
			mdyes160.append(mdmean160[i])
			mdyes098.append(mdmean098[i])
			mdyes105.append(mdmean105[i])
			mdyes110.append(mdmean110[i])
			mdyes140.append(mdmean140[i])
		else:
			mdno125.append(mdmean125[i])
			mdno160.append(mdmean160[i])
			mdno098.append(mdmean098[i])
			mdno105.append(mdmean105[i])
			mdno110.append(mdmean110[i])
			mdno140.append(mdmean140[i])
	mdaryes125 =np.array(mdyes125)
	mdarno125 =np.array(mdno125)
	mdaryes160=np.array(mdyes160)
	mdarno160 =np.array(mdno160)
	mdaryes098=np.array(mdyes098)
	mdarno098 =np.array(mdno098)
	mdaryes105 =np.array(mdyes105)
	mdarno105 =np.array(mdno105)
	mdaryes110=np.array(mdyes110)
	mdarno110 =np.array(mdno110)
	mdaryes140=np.array(mdyes140)
	mdarno140 =np.array(mdno140)
	#change plot opasity to see overlapping point
	#add plot with blob stength/flux(mean) as thierd axis

	mdyesr125_160=mdaryes125/mdaryes160
	#print(len(yesr125_160))
	mdnor125_160=mdarno125/mdarno160
	#print(len(nor125_160))
	mdyesr098_125=mdaryes098/mdaryes125
	mdnor098_125=mdarno098/mdarno125
	plt.rcParams['axes.facecolor'] = 'white'

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(mdyesr125_160,mdyesr098_125, '.', label='Flagged')
	plt.plot(mdnor098_125, mdnor125_160, 's', label='Not Flagged')
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Mean Colors With Flags')	
	plt.legend(prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_colors_flags_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	#27,140,33,25,9,30,148
	blob25 = [mdmean098[26],  mdmean105[26], mdmean110[26], mdmean125[26], mdmean140[26], mdmean160[26]]
	blob100 = [mdmean098[1],mdmean105[1],mdmean110[1],mdmean125[1],mdmean140[1],mdmean160[1]]
	blob65 = [mdmean098[24],  mdmean105[24], mdmean110[24], mdmean125[24], mdmean140[24], mdmean160[24]]
	blob98 = [mdmean098[32],  mdmean105[32], mdmean110[32], mdmean125[32], mdmean140[32], mdmean160[32]]
	blob15 = [mdmean098[8],   mdmean105[8],  mdmean110[8],  mdmean125[8],  mdmean140[8],  mdmean160[8]]
	blob30 = [mdmean098[29],  mdmean105[29], mdmean110[29], mdmean125[29], mdmean140[29], mdmean160[29]]
	blob130 = [mdmean098[14],mdmean105[14],mdmean110[14],mdmean125[14],mdmean140[14],mdmean160[14]]

	lam_list = [9864,10552,11534,12486,13923,15369]




	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(lam_list, blob25, '.', label='blob27')
	plt.plot(lam_list, blob100, '.', label='blob2')
	plt.plot(lam_list, blob65, '.', label='blob25')
	plt.plot(lam_list, blob98, '.', label='blob33')
	plt.plot(lam_list, blob15, '.', label='blob9')
	plt.plot(lam_list, blob30, '.', label='blob30')
	plt.plot(lam_list, blob130, '.', label='blob15')
	plt.xlabel('Pivit Wavelength')
	plt.ylim(0.94,1.02)	
	plt.ylabel('Mean')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Bright_Blob_Mean_sub_set_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	 

	print(len(mdaryes105))
	plt.plot(mdaryes105, mdaryes098 , 'o',label='F098', markersize=4)
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_all_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes105, mdaryes098 , 'o',label='F098', markersize=4)
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_all_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot( mdaryes105, mdaryes160 ,'*',label='F160', markersize=4)
	plt.plot( mdaryes105, mdaryes110 ,'s',label='F110', markersize=4)
	plt.plot( mdaryes105, mdaryes125 ,'d',label='F125', markersize=4)
	plt.plot( mdaryes105, mdaryes140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes110, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes110, mdaryes105 , 's',label='F105', markersize=4)
	plt.plot(mdaryes110, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes110, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_110_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot( mdaryes125, mdaryes160 ,'*',label='F160', markersize=4)
	plt.plot( mdaryes125, mdaryes110 ,'s',label='F110', markersize=4)
	plt.plot( mdaryes125, mdaryes105 ,'d',label='F105', markersize=4)
	plt.plot( mdaryes125, mdaryes140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_125_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes140, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes140, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes140, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes140, mdaryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_140_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes160, mdaryes105 , '*',label='F105', markersize=4)
	plt.plot(mdaryes160, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes160, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes160, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_160_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes110, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes110, mdaryes105 , 's',label='F105', markersize=4)
	plt.plot(mdaryes110, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes110, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_110_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes125, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes125, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes125, mdaryes105 , 'd',label='F105', markersize=4)
	plt.plot(mdaryes125, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_125_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes140, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes140, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes140, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes140, mdaryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_140_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes160, mdaryes105 ,  '*',label='F105', markersize=4)
	plt.plot(mdaryes160, mdaryes110 ,  's',label='F110', markersize=4)
	plt.plot(mdaryes160, mdaryes125 ,  'd',label='F125', markersize=4)
	plt.plot(mdaryes160, mdaryes140 ,  'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_160_md_norm_flaged_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdmean105, mdmean098 , 'o',label='F098', markersize=4)
	plt.plot(mdmean105, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean105, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean105, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean105, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_all_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdmean105, mdmean098 , 'o',label='F098', markersize=4)
	plt.plot(mdmean105, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean105, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean105, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean105, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )	
	plt.ylabel('Mean in each filter')
	plt.title(' Blob Mean per Wavelength')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_all_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot( mdmean105, mdmean160 ,'*',label='F160', markersize=4)
	plt.plot( mdmean105, mdmean110 ,'s',label='F110', markersize=4)
	plt.plot( mdmean105, mdmean125 ,'d',label='F125', markersize=4)
	plt.plot( mdmean105, mdmean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_105_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdmean110, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean110, mdmean105 , 's',label='F105', markersize=4)
	plt.plot(mdmean110, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean110, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_110_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot( mdmean125, mdmean160 ,'*',label='F160', markersize=4)
	plt.plot( mdmean125, mdmean110 ,'s',label='F110', markersize=4)
	plt.plot( mdmean125, mdmean105 ,'d',label='F105', markersize=4)
	plt.plot( mdmean125, mdmean140 ,'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_125_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdmean140, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean140, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean140, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean140, mdmean105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_140_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdmean160, mdmean105 , '*',label='F105', markersize=4)
	plt.plot(mdmean160, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean160, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean160, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.97, 1.0 )
	plt.ylim(0.97, 1.0 )
	plt.ylabel('Mean in each Wide filter')
	plt.title(' Zoomed Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Mean_vs_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot(mdmean105, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean105, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean105, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean105, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 105 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 105')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_105_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdmean110, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean110, mdmean105 , 's',label='F105', markersize=4)
	plt.plot(mdmean110, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean110, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 110 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 110')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_110_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdmean125, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean125, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean125, mdmean105 , 'd',label='F105', markersize=4)
	plt.plot(mdmean125, mdmean140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 125 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 125')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_125_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdmean140, mdmean160 , '*',label='F160', markersize=4)
	plt.plot(mdmean140, mdmean110 , 's',label='F110', markersize=4)
	plt.plot(mdmean140, mdmean125 , 'd',label='F125', markersize=4)
	plt.plot(mdmean140, mdmean105 , 'p',label='F105', markersize=4)
	plt.xlabel('Mean of 140 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 140')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_140_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdmean160, mdmean105 ,  '*',label='F105', markersize=4)
	plt.plot(mdmean160, mdmean110 ,  's',label='F110', markersize=4)
	plt.plot(mdmean160, mdmean125 ,  'd',label='F125', markersize=4)
	plt.plot(mdmean160, mdmean140 ,  'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.ylabel('Mean in each Wide filter')
	plt.title('Blob Mean in each Wide filter Vs. 160')	
	plt.legend(loc = 'best')
	plt.savefig(os.path.join(plot_path,'Blob_Mean_vs_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	#mins
	#filter ratios
	fmin098m = 'F098M_blob_color_min_final_noT0.csv'
	fmin105w = 'F105W_blob_color_min_final_noT0.csv'
	fmin110w = 'F110W_blob_color_min_final_noT0.csv'
	fmin125w = 'F125W_blob_color_min_final_noT0.csv'
	fmin140w = 'F140W_blob_color_min_final_noT0.csv'
	fmin160w = 'F160W_blob_color_min_final_noT0.csv'


	min_098 = pd.read_csv(fmin098m,names=['id', 'min'])
	min_105 = pd.read_csv(fmin105w,names=['id', 'min'])
	min_110 = pd.read_csv(fmin110w,names=['id', 'min'])
	min_125 = pd.read_csv(fmin125w,names=['id', 'min'])
	min_140 = pd.read_csv(fmin140w,names=['id', 'min'])
	min_160 = pd.read_csv(fmin160w,names=['id', 'min'])

	mdminr125_160=(min_125['min']/md125)/(min_160['min']/md160)
	mdminr105_125=(min_105['min']/md105)/(min_125['min']/md125)
	mdminr098_125=(min_098['min']/md098)/(min_125['min']/md125)
	mdminr110_125=(min_110['min']/md110)/(min_125['min']/md125)



	#cm = plt.cm.get_cmap('CMRmap')
	cmap=plt.cm.get_cmap('rainbow', 20)
	plt.rcParams['axes.facecolor'] = 'lightgray'

	z=flux#/(max(flux))#df_160['mean']


	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(mdminr125_160,mdminr105_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('105/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_105_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	sc = plt.scatter(mdminr125_160,mdminr098_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_098_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	sc = plt.scatter(mdminr125_160,mdminr110_125,  c=z ,  cmap=cmap )# label=mf[:5])
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('110/125')
	plt.title(' Blob Min Colors')	
	plt.colorbar(sc)
	plt.clim(0, 180)
	#plt.legend()
	plt.savefig(os.path.join(plot_path,'Blob_Min_colors_110_125v125_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	
	mdmin125 = min_125['min']/md125
	mdmin160 = min_160['min']/md160
	mdmin098 = min_098['min']/md098
	mdmin105 = min_105['min']/md105
	mdmin110 = min_110['min']/md110
	mdmin140 = min_140['min']/md140
	test = min_125['min']
	print('norm =', mdmin125[0], 'ori =', test[0])


	mdyes125 = []
	mdno125 = []
	mdyes160 = []
	mdno160 = []
	mdyes098 = []
	mdno098 = []
	mdyes105 = []
	mdno105 = []
	mdyes110 = []
	mdno110 = []
	mdyes140 = []
	mdno140 = []
	for i in range(len(min125)):
		#print(flags[i])
		if flags[i] == 'Yes':
			mdyes125.append(mdmin125[i])
			mdyes160.append(mdmin160[i])
			mdyes098.append(mdmin098[i])
			mdyes105.append(mdmin105[i])
			mdyes110.append(mdmin110[i])
			mdyes140.append(mdmin140[i])
		else:
			mdno125.append(mdmin125[i])
			mdno160.append(mdmin160[i])
			mdno098.append(mdmin098[i])
			mdno105.append(mdmin105[i])
			mdno110.append(mdmin110[i])
			mdno140.append(mdmin140[i])
	mdaryes125 =np.array(mdyes125)
	mdarno125 =np.array(mdno125)
	mdaryes160=np.array(mdyes160)
	mdarno160 =np.array(mdno160)
	mdaryes098=np.array(mdyes098)
	mdarno098 =np.array(mdno098)
	mdaryes105 =np.array(mdyes105)
	mdarno105 =np.array(mdno105)
	mdaryes110=np.array(mdyes110)
	mdarno110 =np.array(mdno110)
	mdaryes140=np.array(mdyes140)
	mdarno140 =np.array(mdno140)
	#change plot opasity to see overlapping point
	#add plot with blob stength/flux(mean) as thierd axis

	ymdesr125_160=mdaryes125/mdaryes160
	#print(len(yesr125_160))
	mdnor125_160=mdarno125/mdarno160
	#print(len(nor125_160))
	mdyesr098_125=mdaryes098/mdaryes125
	mdnor098_125=mdarno098/mdarno125
	plt.rcParams['axes.facecolor'] = 'white'

	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(mdyesr125_160,mdyesr098_125, '.', label='Flagged')
	plt.plot(mdnor098_125, mdnor125_160, 's', label='Not Flagged')
	plt.xlabel('125/160')
	#plt.xlim(0,0.2)	
	plt.ylabel('098/125')
	plt.title(' Blob Min Colors With Flags')	
	plt.legend()
	plt.savefig(os.path.join(plot_path,'Min_Blob_colors_flags_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	#27,140,33,25,9,30,148
	blob25 = [mdmin098[26],  mdmin105[26],  mdmin110[26], mdmin125[26], mdmin140[26], mdmin160[26]]
	blob100 = [mdmean098[1],mdmean105[1],mdmean110[1],mdmean125[1],mdmean140[1],mdmean160[1]]
	blob65 = [mdmin098[24],  mdmin105[24],  mdmin110[24], mdmin125[24], mdmin140[24], mdmin160[24]]
	blob98 = [mdmin098[32],  mdmin105[32],  mdmin110[32], mdmin125[32], mdmin140[32], mdmin160[32]]
	blob15 = [mdmin098[8],   mdmin105[8],   mdmin110[8],  mdmin125[8],  mdmin140[8],  mdmin160[8]]
	blob30 = [mdmin098[29],  mdmin105[29],  mdmin110[29], mdmin125[29], mdmin140[29], mdmin160[29]]
	blob130 = [mean098[14],mean105[14],mean110[14],mean125[14],mean140[14],mean160[14]]

	lam_list = [9864,10552,11534,12486,13923,15369]




	#plt.plot(df_mean['mean'],df_sum['sum'], '.', label=mf[:5])
	plt.plot(lam_list, blob25, '.', label='blob27')
	plt.plot(lam_list, blob100, '.', label='blob140')
	plt.plot(lam_list, blob65, '.', label='blob25')
	plt.plot(lam_list, blob98, '.', label='blob33')
	plt.plot(lam_list, blob15, '.', label='blob9')
	plt.plot(lam_list, blob30, '.', label='blob30')
	plt.plot(lam_list, blob130, '.', label='blob148')
	plt.xlabel('Pivit Wavelength')
	#plt.xlim(0,0.2)	
	plt.ylabel('Min')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Bright_Blob_Min_sub_set_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	 

	print(len(mdaryes105))
	plt.plot(mdaryes105, mdaryes098 , 'o',label='F098', markersize=4)
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	#plt.xlim(0,0.2)	
	plt.ylabel('Min in each filter')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_105_all_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes105, mdaryes098, 'o',label='F098', markersize=4)
	plt.plot(mdaryes105, mdaryes160, '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110, 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125, 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140, 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )	
	plt.ylabel('Min in each filter')
	plt.title(' Blob Min per Wavelength')	
	plt.legend(loc = 'best', prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_105_all_md_norm_noT0_v4_no_norm.png'))
	plt.clf()



	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	print(len(mdaryes160),len(mdaryes105),len(mdaryes110),len(mdaryes125),len(mdaryes140))
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 105')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_105_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes110, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes110, mdaryes105 , 's',label='F105', markersize=4)
	plt.plot(mdaryes110, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes110, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 110 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 110')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_110_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes125, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes125, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes125, mdaryes105 , 'd',label='F105', markersize=4)
	plt.plot(mdaryes125, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 125 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 125')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_125_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes140, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes140, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes140, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes140, mdaryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Min of 140 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 140')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_140_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes160, mdaryes105 , '*',label='F105', markersize=4)
	plt.plot(mdaryes160, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes160, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes160, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Mean of 160 filter')
	plt.xlim(0.85, 1.0 )
	plt.ylim(0.85, 1.0 )
	plt.ylabel('Min in each Wide filter')
	plt.title(' Zoomed Blob Min in each Wide filter Vs. 160')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Zoomed_Blob_Min_vs_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	#plt.plot(aryes098 , aryes105, 'o',label='F098')
	plt.plot(mdaryes105, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes105, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes105, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes105, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 105 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 105')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_105_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes110, mdaryes160, '*',label='F160', markersize=4)
	plt.plot(mdaryes110, mdaryes105, 's',label='F105', markersize=4)
	plt.plot(mdaryes110, mdaryes125, 'd',label='F125', markersize=4)
	plt.plot(mdaryes110, mdaryes140, 'p',label='F140', markersize=4)
	plt.xlabel('Min of 110 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 110')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_110_md_norm_noT0_v4_no_norm.png'))
	plt.clf()

	plt.plot(mdaryes125, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes125, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes125, mdaryes105 , 'd',label='F105', markersize=4)
	plt.plot(mdaryes125, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 125 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 125')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_125_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes140, mdaryes160 , '*',label='F160', markersize=4)
	plt.plot(mdaryes140, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes140, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes140, mdaryes105 , 'p',label='F105', markersize=4)
	plt.xlabel('Min of 140 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 140')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_140_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


	plt.plot(mdaryes160, mdaryes105 , '*',label='F105', markersize=4)
	plt.plot(mdaryes160, mdaryes110 , 's',label='F110', markersize=4)
	plt.plot(mdaryes160, mdaryes125 , 'd',label='F125', markersize=4)
	plt.plot(mdaryes160, mdaryes140 , 'p',label='F140', markersize=4)
	plt.xlabel('Min of 160 filter')
	plt.ylabel('Min in each Wide filter')
	plt.title('Blob Min in each Wide filter Vs. 160')	
	plt.legend(loc = 'best',prop={'size': 6})
	plt.savefig(os.path.join(plot_path,'Blob_Min_vs_160_md_norm_noT0_v4_no_norm.png'))
	plt.clf()


main()