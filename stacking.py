import numpy as np
import matplotlib
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import matplotlib.patches as mpatches
import os
from multiprocessing import Pool
import argparse
import sys



#-------------------------
# Median Stack
#-------------------------
def median_stack1(j):
	"""This is a function to create a median stack of a data cube through the z axis. It will print the column it is currently working on.

	median_stack1(j)
	Parameters:
		j: int
			The number of columns to median stack together at a time


	"""
	col_median1=np.nanmedian(data_cube1[:,:,j],axis=0)
#	print(j)
	return col_median1
	
def median_stack2(j):
	"""This is a function to create a median stack of a data cube through the z axis. It will print the column it is currently working on.

	median_stack1(j)
	Parameters:
		j: int
			The number of columns to median stack together at a time


	"""
	col_median2=np.nanmedian(data_cube2[:,:,j],axis=0)
#	print(j)
	return col_median2
#-------------------------
# Mean Stack
#-------------------------	
def mean_stack1(j):
	"""This is a function to create a mean stack of a data cube through the z axis.

	mean_stack1(j)
	Parameters:
		j: int
			The number of columns to mean stack together at a time


	"""
	col_mean1=np.nanmean(data_cube1[:,:,j],axis=0)
	#print(j)
	return col_mean1


def mean_stack2(j):
	"""This is a function to create a mean stack of a data cube through the z axis.

	mean_stack1(j)
	Parameters:
		j: int
			The number of columns to mean stack together at a time


	"""
	col_mean2=np.nanmean(data_cube2[:,:,j],axis=0)
	#print(j)
	return col_mean2

def sum_stacked1(j):
	"""This is a function to create a sum stack of a data cube through the z axis.

	sum_stacked1(j)
	Parameters:
		j: int
			The number of columns to sum stack together at a time


	"""
	col_err1=np.nansum(err_data_cube1[:,:,j]**2,axis=0)
#	nf1=np.count_nonzero(~np.isnan(err_data_cube1[:,:,j]))
	#print(j)
	return col_err1


def sum_stacked2(j):
	"""This is a function to create a mean stack of a data cube through the z axis.

	mean_stack1(j)
	Parameters:
		j: int
			The number of columns to mean stack together at a time


	"""
	col_err2=np.nansum(err_data_cube2[:,:,j]**2,axis=0)
	#print(j)
	return col_err2		


def stacking(path):
    """This function averages all science extentsions and propagates the error extentsions 

    Author: Myles McKay
    Date: April 30, 2017

    Function:
        stacking(path)

    Parameters:
        path: str
            The path to the directory containing raw.fits files.

    Output:
        Averaged bias image and propagated error arrays
    Example:

        python stacking.py --path='/bias/flt/file/location/' 

    """

	#Masked the data with the DQ array 
	i=0
	for fitsname in files:
	    h = fits.open(fitsname)
	    sci_chip1=h[4].data
	    err_chip1=h[5].data
	    dq_chip1=h[6].data
	    sci_chip2=h[1].data
	    err_chip2=h[2].data
	    dq_chip2=h[3].data
	    h.close()

	    sci_chip1[dq_chip1 !=0]=np.nan
	    sci_chip2[dq_chip2 !=0]=np.nan
	    err_chip1[dq_chip1 !=0]=np.nan
	    err_chip2[dq_chip2 !=0]=np.nan
	
	
		#Inputs the masked data into 3D array
	    data_cube1[i] = sci_chip1
	    data_cube2[i] = sci_chip2
	    err_data_cube1[i]=err_chip1
	    err_data_cube2[i]=err_chip2
	    i+=1
	    print(i)
	
	
	nf1=np.count_nonzero(~np.isnan(err_data_cube1[:,:,range(4096)]),axis=0)
	nf2=np.count_nonzero(~np.isnan(err_data_cube2[:,:,range(4096)]),axis=0)
	#Using parellel computing to median stack the columns for faster results (Do not use 8 when running the locally)	
	p=Pool(20)
	
	#result1 = p.map(median_stack1,range(4096))
	#result2 = p.map(median_stack2,range(4096))
	result1 = p.map(mean_stack1,range(4096))
	result2 = p.map(mean_stack2,range(4096))
	err_result1=p.map(sum_stacked1,range(4096))
	err_result2=p.map(sum_stacked2,range(4096))
	
	
	#Puts the final 2d list into a numpy array
	result1=np.array(result1)
	result2=np.array(result2)
	err_result1=np.array(err_result1)
	err_result2=np.array(err_result2)
	
	#print(err_result1[1344:1350,1303:1309])
	
	#print(err_result1[500,500])
	#print(err_result2[500,500])
	#The data deminsions changes when the median stack occures so we transpose the data to re shape the array
	result1 =np.transpose(result1)
	result2 =np.transpose(result2)
	err_result1=np.transpose(err_result1)
	err_result2=np.transpose(err_result2)
	#print(err_result2[1703:1708,2693:2698])
	#nf1=np.count_nonzero(~np.isnan(err_result1))
	#nf2=np.count_nonzero(~np.isnan(err_result2))
	print(nf1[1703:1708,2693:2698])
	print(nf2[1703:1708,2693:2698])
	err_result1=(1/nf1)*np.sqrt(err_result1)
	err_result2=(1/nf2)*np.sqrt(err_result2)
	
	#print(err_result1[1344:1350,1303:1309])
	print(err_result2[1703:1708,2693:2698])
	#print(err_result1[500,500])
	#print(err_result2[500,500])
	#print (result[0:3,0:3])
	#prints the shape of the array to confirm the disired deminsions
	#print (np.shape(result1))
	#print (np.shape(result2))
	#print(result1[0,0])
	#print(result2[0,0])
	#Writes a new fits for the data
	new_hdul = fits.HDUList()
	new_hdul.append(fits.ImageHDU(result1))
	new_hdul.append(fits.ImageHDU(err_result1))
	new_hdul.append(fits.ImageHDU(result2))
	new_hdul.append(fits.ImageHDU(err_result2))
	new_hdul.writeto('bias_crr_err_stacked_files.fits',overwrite=True)
	
	#plt.savefig('median_stack_bias1.pdf')
	#plt.show()
	
	#moves new fits files to a new directory
	#os.system('mv Master_flc_mean_stack_chip*.fits /grp/hst/wfc3v/mmckay/2009biasfiles/median_stack')



def parse_args():
    """Parses command line arguments.

    Parameters:
        nothing

    Returns:
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    Outputs:
        nothing
    """

    #Create help string:
    path_help = 'Path to the folder with files to run tweakreg.'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-path', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)


    # Parse args:
    args = parser.parse_args()

    return args
# -------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    os.chdir(args.path)
    files= glob.glob('*flt.fits')
    hdr = fits.getheader(files[0], 1)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nf=len(files)
    data_cube1 = np.zeros((nf, ny, nx), dtype=float)
    data_cube2 = np.zeros((nf, ny, nx), dtype=float)
    err_data_cube1 = np.zeros((nf, ny, nx), dtype=float)
    err_data_cube2 = np.zeros((nf, ny, nx), dtype=float)
    stacking(args.path)













   