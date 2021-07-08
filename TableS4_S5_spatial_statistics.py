#!/usr/bin/env python2.7

import sys
import os, glob
import numpy as np 
import pandas as pd

def calc_spatial_sig(in_file):
	data = pd.read_csv(in_file,index_col=None,header=None, sep=',',parse_dates=True, squeeze=True)                          
	lower = data.values

	# get domain mask
	oceanmask = np.genfromtxt('oceanmask.csv',delimiter=',')
	oceanmask[0:23,:] = np.nan
	oceanmask[98:116,:] = np.nan
	oceanmask[:,100:117] = np.nan

	sig = lower
	sig[sig<0.001]=0
	sig[sig>0.000999999]=1
	sig = sig+oceanmask

	sig_ratio = np.nansum(sig)/3565

	return sig_ratio

def calc_region_mean(in_file1,in_file2,ptile):
	# in_file1 : natural forcing simulations
	# in_file2 : full forcing simulations
	# ptile    : percentile 

	# get ocean mask
	oceanmask = np.genfromtxt('oceanmask.csv',delimiter=',')
	oceanmask[0:23,:] = np.nan
	oceanmask[98:116,:] = np.nan
	oceanmask[:,100:117] = np.nan

    # readin data
	data = pd.read_csv(in_file1,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	natdata = data.values
	natptile = np.percentile(natdata,ptile,axis=1)
	natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

	data = pd.read_csv(in_file2,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	fuldata = data.values
	fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

	delta = np.divide(fulcount,natcount)-1
	fwi_domain = np.reshape(delta, (116,110))
	fwi_domain = fwi_domain+oceanmask

	region_mean = np.nanmean(fwi_domain)

	return region_mean

###############################

out_file = 'Table_S4.csv'
outdata = np.empty([9,2])

############ fwi  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'FWI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'fwi_domain_confint_lower.txt')

outdata[0,0] = calc_spatial_sig(in_file3)
outdata[0,1] = calc_region_mean(in_file1,in_file2,ptile)

############ isi  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'ISI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'isi_domain_confint_lower.txt')

outdata[1,0] = calc_spatial_sig(in_file3)
outdata[1,1]= calc_region_mean(in_file1,in_file2,ptile)

############ bui  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'BUI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'bui_domain_confint_lower.txt')

outdata[2,0] = calc_spatial_sig(in_file3)
outdata[2,1] = calc_region_mean(in_file1,in_file2,ptile)

############ hdw  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'HDW_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'hdw_domain_confint_lower.txt')

outdata[3,0] = calc_spatial_sig(in_file3)
outdata[3,1] = calc_region_mean(in_file1,in_file2,ptile)

############ fosberg  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'FOS_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'fos_domain_confint_lower.txt')

outdata[4,0] = calc_spatial_sig(in_file3)
outdata[4,1] = calc_region_mean(in_file1,in_file2,ptile)

############ tmax  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'TMAX_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'tmax_domain_confint_lower.txt')

outdata[5,0] = calc_spatial_sig(in_file3)
outdata[5,1] = calc_region_mean(in_file1,in_file2,ptile)

############ rh  ##########
ptile = 5.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'RH_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'rh_domain_confint_lower.txt')

outdata[6,0] = calc_spatial_sig(in_file3)
outdata[6,1] = calc_region_mean(in_file1,in_file2,ptile)

############ vpd  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'VPD_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'vpd_domain_confint_lower.txt')

outdata[7,0] = calc_spatial_sig(in_file3)
outdata[7,1] = calc_region_mean(in_file1,in_file2,ptile)

############ map windspeed  ##########
ptile = 95.0
in_dir = './data/'
in_file1 = os.path.join(in_dir+'WIND_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'wind_domain_confint_lower.txt')

outdata[8,0] = calc_spatial_sig(in_file3)
outdata[8,1] = calc_region_mean(in_file1,in_file2,ptile)

#####################################
print outdata
np.savetxt(out_file,outdata,fmt='%6.4f')

############################################
############### table S5 ###################
############################################

def get_mask(in_file,ptile):
	df = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)
	datptile = np.percentile(df.values,ptile,axis=1)
	datptile_mat = np.transpose([datptile]*2000)
	mask = df[df>datptile_mat]
	mask[mask>-10000] = 1

	return mask

################# FWI mask #################
ptile = 95.0
in_dir = './data/'
in_file = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
fwi95_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'FWI_allsets_b840_domain_byFWI.txt')
fwi95nat_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
isi95_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'ISI_allsets_b840_domain_byFWI.txt')
isi95nat_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
bui95_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'BUI_allsets_b840_domain_byFWI.txt')
bui95nat_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
hdw95_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'HDW_allsets_b840_domain_byFWI.txt')
hdw95nat_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
ffwi95_mask = get_mask(in_file,ptile)

in_file = os.path.join(in_dir+'FOS_allsets_b840_domain_byFWI.txt')
ffwi95nat_mask = get_mask(in_file,ptile)


########## tmax ############

var_diff = np.empty([5,4])
in_dir = './data/'
in_file = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
var = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

in_file = os.path.join(in_dir+'TMAX_allsets_b840_domain_byFWI.txt')
var_nat = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

var95 = var*fwi95_mask
var95nat = var_nat*fwi95nat_mask
var_diff[0,0] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*isi95_mask
var95nat = var_nat*isi95nat_mask
var_diff[1,0] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*bui95_mask
var95nat = var_nat*bui95nat_mask
var_diff[2,0] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*hdw95_mask
var95nat = var_nat*hdw95nat_mask
var_diff[3,0] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*ffwi95_mask
var95nat = var_nat*ffwi95nat_mask
var_diff[4,0] = np.nanmean(var95) - np.nanmean(var95nat)

############## RH ###############

in_dir = './data/'
in_file = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
var = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

in_file = os.path.join(in_dir+'RH_allsets_b840_domain_byFWI.txt')
var_nat = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

var95 = var*fwi95_mask
var95nat = var_nat*fwi95nat_mask
var_diff[0,1] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*isi95_mask
var95nat = var_nat*isi95nat_mask
var_diff[1,1] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*bui95_mask
var95nat = var_nat*bui95nat_mask
var_diff[2,1] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*hdw95_mask
var95nat = var_nat*hdw95nat_mask
var_diff[3,1] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*ffwi95_mask
var95nat = var_nat*ffwi95nat_mask
var_diff[4,1] = np.nanmean(var95) - np.nanmean(var95nat)

############## vpd ##############

in_dir = './data/'
in_file = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
var = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

in_file = os.path.join(in_dir+'VPD_allsets_b840_domain_byFWI.txt')
var_nat = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

var95 = var*fwi95_mask
var95nat = var_nat*fwi95nat_mask
var_diff[0,2] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*isi95_mask
var95nat = var_nat*isi95nat_mask
var_diff[1,2] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*bui95_mask
var95nat = var_nat*bui95nat_mask
var_diff[2,2] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*hdw95_mask
var95nat = var_nat*hdw95nat_mask
var_diff[3,2] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*ffwi95_mask
var95nat = var_nat*ffwi95nat_mask
var_diff[4,2] = np.nanmean(var95) - np.nanmean(var95nat)

########## windspeed ############

in_dir = './data/'
in_file = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
var = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

in_file = os.path.join(in_dir+'WIND_allsets_b840_domain_byFWI.txt')
var_nat = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)

var95 = var*fwi95_mask
var95nat = var_nat*fwi95nat_mask
var_diff[0,3] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*isi95_mask
var95nat = var_nat*isi95nat_mask
var_diff[1,3] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*bui95_mask
var95nat = var_nat*bui95nat_mask
var_diff[2,3] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*hdw95_mask
var95nat = var_nat*hdw95nat_mask
var_diff[3,3] = np.nanmean(var95) - np.nanmean(var95nat)

var95 = var*ffwi95_mask
var95nat = var_nat*ffwi95nat_mask
var_diff[4,3] = np.nanmean(var95) - np.nanmean(var95nat)

np.savetxt('Table_S5.csv',var_diff,fmt='%6.4f')
