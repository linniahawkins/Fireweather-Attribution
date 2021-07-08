#!/usr/bin/env python2.7

import sys
import os, glob
from netCDF4 import Dataset
import numpy as np 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
from cartopy import config
import cartopy.feature as cfeature
from scipy import stats

# readin domain masks

oregon_mask = np.genfromtxt('PSA_oregon_mask.csv',delimiter=',')
oregon_mask[oregon_mask == 0] = np.nan 
tmp2 = np.reshape(oregon_mask,12760)
oregon_cells = np.argwhere(tmp2>0).flatten()

paradise_mask = np.genfromtxt('PSA_CAnorthsierra_mask.csv',delimiter=',')
paradise_mask[paradise_mask == 0] = np.nan 
tmp2 = np.reshape(paradise_mask,12760)
paradise_cells = np.argwhere(tmp2>0).flatten()

tubbs_mask = np.genfromtxt('PSA_CAnorthcoast_mask.csv',delimiter=',')
tubbs_mask[tubbs_mask == 0] = np.nan 
tmp2 = np.reshape(tubbs_mask,12760)
tubbs_cells = np.argwhere(tmp2>0).flatten()

malibu_mask = np.genfromtxt('PSA_CAsouthcoast_mask.csv',delimiter=',')
malibu_mask[malibu_mask == 0] = np.nan 
tmp2 = np.reshape(malibu_mask,12760)
malibu_cells = np.argwhere(tmp2>0).flatten()


ptile = 95.0

def bootstrap_confint(em1,em2, ptile, direction="above", c=[5, 95], bsn=1e5):
    # ptile = percentile 
    # direction = "above" find frequency of occurance above ptile (or below)
    # c = confidence intervals (percentiles) to calculate
    # bsn = boot strap number, number of times to resample the distribution
    nat_data = em1 # 9x2000
    ful_data = em2
    # create the store
    sample_store = np.zeros((nat_data.shape[0],int(bsn)), 'f')
    # define sample size
    ssize = nat_data.shape[1]
    # do the resampling
    for s in range(0, int(bsn)):
        x1 = np.random.uniform(0,nat_data.shape[1],ssize).astype(int)
        x2 = np.random.uniform(0,nat_data.shape[1],ssize).astype(int)
        n_data = nat_data[:,x1]
        f_data = ful_data[:,x2]
        val = np.percentile(n_data,ptile,axis=1)
        nat_occ = []; ful_occ = [];
        if direction == "above":
            for ij in range(len(val)):
                nat_occ.append(float(np.count_nonzero(n_data[ij,:] > val[ij])))
                ful_occ.append(float(np.count_nonzero(f_data[ij,:] > val[ij])))
        elif direction == "below":
            for ij in range(len(val)):
                nat_occ.append(float(np.count_nonzero(n_data[ij,:] < val[ij])))
                ful_occ.append(float(np.count_nonzero(f_data[ij,:] < val[ij])))
        else:
            print 'check direction input'
        sample_store[:,s] = 100*(np.divide(ful_occ,nat_occ)-1)

    # now for each confidence interval find the frequency above/below ptile
    spatial_average = np.nanmean(sample_store,0)
    low = np.percentile(spatial_average, c[0],0)
    high = np.percentile(spatial_average, c[1],0)
    conf_inter = [low,high]
    
    return conf_inter

####################################
# readin fwi data
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'FWI_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
FWI_tubbs_nat = natdata[tubbs_cells,:]
FWI_paradise_nat = natdata[paradise_cells,:]
FWI_oregon_nat = natdata[oregon_cells,:]
FWI_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
FWI_tubbs_ful = fuldata[tubbs_cells,:]
FWI_paradise_ful = fuldata[paradise_cells,:]
FWI_oregon_ful = fuldata[oregon_cells,:]
FWI_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
fwi_domain = np.reshape(delta, (116,110))
FWI_tubbs = np.nanmean(fwi_domain*tubbs_mask)
FWI_paradise = np.nanmean(fwi_domain*paradise_mask)
FWI_oregon = np.nanmean(fwi_domain*oregon_mask)
FWI_malibu = np.nanmean(fwi_domain*malibu_mask)

fwi_tubbs_confint = bootstrap_confint(FWI_tubbs_nat,FWI_tubbs_ful,ptile,direction="above")
fwi_paradise_confint = bootstrap_confint(FWI_paradise_nat,FWI_paradise_ful,ptile,direction="above")
fwi_oregon_confint = bootstrap_confint(FWI_oregon_nat,FWI_oregon_ful,ptile,direction="above")
fwi_malibu_confint = bootstrap_confint(FWI_malibu_nat,FWI_malibu_ful,ptile,direction="above")

####################################
# readin isi data
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'ISI_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
ISI_tubbs_nat = natdata[tubbs_cells,:]
ISI_paradise_nat = natdata[paradise_cells,:]
ISI_oregon_nat = natdata[oregon_cells,:]
ISI_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
ISI_tubbs_ful = fuldata[tubbs_cells,:]
ISI_paradise_ful = fuldata[paradise_cells,:]
ISI_oregon_ful = fuldata[oregon_cells,:]
ISI_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
isi_domain = np.reshape(delta, (116,110))
ISI_tubbs = np.nanmean(isi_domain*tubbs_mask)
ISI_paradise = np.nanmean(isi_domain*paradise_mask)
ISI_oregon = np.nanmean(isi_domain*oregon_mask)
ISI_malibu = np.nanmean(isi_domain*malibu_mask)

isi_tubbs_confint = bootstrap_confint(ISI_tubbs_nat,ISI_tubbs_ful,ptile,direction="above")
isi_paradise_confint = bootstrap_confint(ISI_paradise_nat,ISI_paradise_ful,ptile,direction="above")
isi_oregon_confint = bootstrap_confint(ISI_oregon_nat,ISI_oregon_ful,ptile,direction="above")
isi_malibu_confint = bootstrap_confint(ISI_malibu_nat,ISI_malibu_ful,ptile,direction="above")

####################################
# readin bui counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'BUI_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
BUI_tubbs_nat = natdata[tubbs_cells,:]
BUI_paradise_nat = natdata[paradise_cells,:]
BUI_oregon_nat = natdata[oregon_cells,:]
BUI_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
BUI_tubbs_ful = fuldata[tubbs_cells,:]
BUI_paradise_ful = fuldata[paradise_cells,:]
BUI_oregon_ful = fuldata[oregon_cells,:]
BUI_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
bui_domain = np.reshape(delta, (116,110))
BUI_tubbs = np.nanmean(bui_domain*tubbs_mask)
BUI_paradise = np.nanmean(bui_domain*paradise_mask)
BUI_oregon = np.nanmean(bui_domain*oregon_mask)
BUI_malibu = np.nanmean(bui_domain*malibu_mask)

bui_tubbs_confint = bootstrap_confint(BUI_tubbs_nat,BUI_tubbs_ful,ptile,direction="above")
bui_paradise_confint = bootstrap_confint(BUI_paradise_nat,BUI_paradise_ful,ptile,direction="above")
bui_oregon_confint = bootstrap_confint(BUI_oregon_nat,BUI_oregon_ful,ptile,direction="above")
bui_malibu_confint = bootstrap_confint(BUI_malibu_nat,BUI_malibu_ful,ptile,direction="above")

####################################
# readin HDW counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'HDW_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
HDW_tubbs_nat = natdata[tubbs_cells,:]
HDW_paradise_nat = natdata[paradise_cells,:]
HDW_oregon_nat = natdata[oregon_cells,:]
HDW_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
HDW_tubbs_ful = fuldata[tubbs_cells,:]
HDW_paradise_ful = fuldata[paradise_cells,:]
HDW_oregon_ful = fuldata[oregon_cells,:]
HDW_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
hdw_domain = np.reshape(delta, (116,110))
HDW_tubbs = np.nanmean(hdw_domain*tubbs_mask)
HDW_paradise = np.nanmean(hdw_domain*paradise_mask)
HDW_oregon = np.nanmean(hdw_domain*oregon_mask)
HDW_malibu = np.nanmean(hdw_domain*malibu_mask)

hdw_tubbs_confint = bootstrap_confint(HDW_tubbs_nat,HDW_tubbs_ful,ptile,direction="above")
hdw_paradise_confint = bootstrap_confint(HDW_paradise_nat,HDW_paradise_ful,ptile,direction="above")
hdw_oregon_confint = bootstrap_confint(HDW_oregon_nat,HDW_oregon_ful,ptile,direction="above")
hdw_malibu_confint = bootstrap_confint(HDW_malibu_nat,HDW_malibu_ful,ptile,direction="above")

####################################
# readin Fosberg counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'FOS_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
FOS_tubbs_nat = natdata[tubbs_cells,:]
FOS_paradise_nat = natdata[paradise_cells,:]
FOS_oregon_nat = natdata[oregon_cells,:]
FOS_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]


in_file = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
FOS_tubbs_ful = fuldata[tubbs_cells,:]
FOS_paradise_ful = fuldata[paradise_cells,:]
FOS_oregon_ful = fuldata[oregon_cells,:]
FOS_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
fos_domain = np.reshape(delta, (116,110))
FOS_tubbs = np.nanmean(fos_domain*tubbs_mask)
FOS_paradise = np.nanmean(fos_domain*paradise_mask)
FOS_oregon = np.nanmean(fos_domain*oregon_mask)
FOS_malibu = np.nanmean(fos_domain*malibu_mask)

fos_tubbs_confint = bootstrap_confint(FOS_tubbs_nat,FOS_tubbs_ful,ptile,direction="above")
fos_paradise_confint = bootstrap_confint(FOS_paradise_nat,FOS_paradise_ful,ptile,direction="above")
fos_oregon_confint = bootstrap_confint(FOS_oregon_nat,FOS_oregon_ful,ptile,direction="above")
fos_malibu_confint = bootstrap_confint(FOS_malibu_nat,FOS_malibu_ful,ptile,direction="above")

####################################
# readin Tmax counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'TMAX_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
TMAX_tubbs_nat = natdata[tubbs_cells,:]
TMAX_paradise_nat = natdata[paradise_cells,:]
TMAX_oregon_nat = natdata[oregon_cells,:]
TMAX_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
TMAX_tubbs_ful = fuldata[tubbs_cells,:]
TMAX_paradise_ful = fuldata[paradise_cells,:]
TMAX_oregon_ful = fuldata[oregon_cells,:]
TMAX_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
tmax_domain = np.reshape(delta, (116,110))
TMAX_tubbs = np.nanmean(tmax_domain*tubbs_mask)
TMAX_paradise = np.nanmean(tmax_domain*paradise_mask)
TMAX_oregon = np.nanmean(tmax_domain*oregon_mask)
TMAX_malibu = np.nanmean(tmax_domain*malibu_mask)

tmax_tubbs_confint = bootstrap_confint(TMAX_tubbs_nat,TMAX_tubbs_ful,ptile,direction="above")
tmax_paradise_confint = bootstrap_confint(TMAX_paradise_nat,TMAX_paradise_ful,ptile,direction="above")
tmax_oregon_confint = bootstrap_confint(TMAX_oregon_nat,TMAX_oregon_ful,ptile,direction="above")
tmax_malibu_confint = bootstrap_confint(TMAX_malibu_nat,TMAX_malibu_ful,ptile,direction="above")

####################################
# readin RH counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'RH_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
RH_tubbs_nat = natdata[tubbs_cells,:]
RH_paradise_nat = natdata[paradise_cells,:]
RH_oregon_nat = natdata[oregon_cells,:]
RH_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,100-ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] < natptile[cell])) for cell in range(natdata.shape[0])]


in_file = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
RH_tubbs_ful = fuldata[tubbs_cells,:]
RH_paradise_ful = fuldata[paradise_cells,:]
RH_oregon_ful = fuldata[oregon_cells,:]
RH_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] < natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
rh_domain = np.reshape(delta, (116,110))
RH_tubbs = np.nanmean(rh_domain*tubbs_mask)
RH_paradise = np.nanmean(rh_domain*paradise_mask)
RH_oregon = np.nanmean(rh_domain*oregon_mask)
RH_malibu = np.nanmean(rh_domain*malibu_mask)

rh_tubbs_confint = bootstrap_confint(RH_tubbs_nat,RH_tubbs_ful,100-ptile,direction="below")
rh_paradise_confint = bootstrap_confint(RH_paradise_nat,RH_paradise_ful,100-ptile,direction="below")
rh_oregon_confint = bootstrap_confint(RH_oregon_nat,RH_oregon_ful,100-ptile,direction="below")
rh_malibu_confint = bootstrap_confint(RH_malibu_nat,RH_malibu_ful,100-ptile,direction="below")


####################################
# readin VPD counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'VPD_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
VPD_tubbs_nat = natdata[tubbs_cells,:]
VPD_paradise_nat = natdata[paradise_cells,:]
VPD_oregon_nat = natdata[oregon_cells,:]
VPD_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]


in_file = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
VPD_tubbs_ful = fuldata[tubbs_cells,:]
VPD_paradise_ful = fuldata[paradise_cells,:]
VPD_oregon_ful = fuldata[oregon_cells,:]
VPD_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
vpd_domain = np.reshape(delta, (116,110))
VPD_tubbs = np.nanmean(vpd_domain*tubbs_mask)
VPD_paradise = np.nanmean(vpd_domain*paradise_mask)
VPD_oregon = np.nanmean(vpd_domain*oregon_mask)
VPD_malibu = np.nanmean(vpd_domain*malibu_mask)

vpd_tubbs_confint = bootstrap_confint(VPD_tubbs_nat,VPD_tubbs_ful,ptile,direction="above")
vpd_paradise_confint = bootstrap_confint(VPD_paradise_nat,VPD_paradise_ful,ptile,direction="above")
vpd_oregon_confint = bootstrap_confint(VPD_oregon_nat,VPD_oregon_ful,ptile,direction="above")
vpd_malibu_confint = bootstrap_confint(VPD_malibu_nat,VPD_malibu_ful,ptile,direction="above")

####################################
# readin wind counts
####################################
in_dir = '/Users/linniahawkins/Desktop/campfire/data/final/'
in_file = os.path.join(in_dir+'WIND_allsets_b840_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
natdata = data.values
WS_tubbs_nat = natdata[tubbs_cells,:]
WS_paradise_nat = natdata[paradise_cells,:]
WS_oregon_nat = natdata[oregon_cells,:]
WS_malibu_nat = natdata[malibu_cells,:]
natptile = np.percentile(natdata,ptile,axis=1)
natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

in_file = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
data = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
fuldata = data.values
WS_tubbs_ful = fuldata[tubbs_cells,:]
WS_paradise_ful = fuldata[paradise_cells,:]
WS_oregon_ful = fuldata[oregon_cells,:]
WS_malibu_ful = fuldata[malibu_cells,:]

fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]
delta = 100*(np.divide(fulcount,natcount)-1)
ws_domain = np.reshape(delta, (116,110))
WS_tubbs = np.nanmean(ws_domain*tubbs_mask)
WS_paradise = np.nanmean(ws_domain*paradise_mask)
WS_oregon = np.nanmean(ws_domain*oregon_mask)
WS_malibu = np.nanmean(ws_domain*malibu_mask)

ws_tubbs_confint = bootstrap_confint(WS_tubbs_nat,WS_tubbs_ful,ptile,direction="above")
ws_paradise_confint = bootstrap_confint(WS_paradise_nat,WS_paradise_ful,ptile,direction="above")
ws_oregon_confint = bootstrap_confint(WS_oregon_nat,WS_oregon_ful,ptile,direction="above")
ws_malibu_confint = bootstrap_confint(WS_malibu_nat,WS_malibu_ful,ptile,direction="above")

####### Bar chart Tubbs ###########################
plt.style.use('seaborn-deep')

lower_idx = [e2-e1 for (e1, e2) in zip([fwi_tubbs_confint[0],isi_tubbs_confint[0],bui_tubbs_confint[0],hdw_tubbs_confint[0],fos_tubbs_confint[0]],[FWI_tubbs,ISI_tubbs,BUI_tubbs,HDW_tubbs,FOS_tubbs])]
upper_idx = [e1-e2 for (e1, e2) in zip([fwi_tubbs_confint[1],isi_tubbs_confint[1],bui_tubbs_confint[1],hdw_tubbs_confint[1],fos_tubbs_confint[1]],[FWI_tubbs,ISI_tubbs,BUI_tubbs,HDW_tubbs,FOS_tubbs])]

lower_met = [e2-e1 for (e1,e2) in zip([tmax_tubbs_confint[0],vpd_tubbs_confint[0],rh_tubbs_confint[0],ws_tubbs_confint[0]],[TMAX_tubbs,VPD_tubbs,RH_tubbs,WS_tubbs])]
upper_met = [e1-e2 for (e1,e2) in zip([tmax_tubbs_confint[1],vpd_tubbs_confint[1],rh_tubbs_confint[1],ws_tubbs_confint[1]],[TMAX_tubbs,VPD_tubbs,RH_tubbs,WS_tubbs])]

plt.figure(num=None, figsize=(17, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.9, wspace=0.3, hspace=0.3)
plt.rcParams.update({'font.size': 20})
ax1 = plt.subplot(1,1,1)
ax1.bar([6,7,8,9], [TMAX_tubbs,VPD_tubbs,RH_tubbs,WS_tubbs],yerr=[lower_met,upper_met], width=0.8,color='darkgrey')
ax1.bar([1,2,3,4,5], [FWI_tubbs,ISI_tubbs,BUI_tubbs,HDW_tubbs,FOS_tubbs],yerr=[lower_idx,upper_idx],width=0.8,color='darkgrey')
ax1.set_ylim([-20,130])
ax1.set_ylabel('change in frequency')
ax2= ax1.twinx()
ax2.plot([0,10],[0,0],color='k')
#ax2.plot([0,10],[0.05,0.05],color='k',linestyle='--')
ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel('return period')
ax2.set_yticks([0.0,0.3333,1])
ax2.set_yticklabels(['1/20','1/15','1/10'])

plt.xticks([1,2,3,4,5,6,7,8,9],['FWI','ISI','BUI','HDW','FFWI','TA','VPD','RH','WS'])

plt.title('Tubbs Fire',fontsize=24)
plt.box(False)

plotname = 'Tubbs_barchart_PSAmean_deltaFrequency'
plt.savefig(plotname, dpi=300)

####### Bar chart Paradise ###########################
plt.style.use('seaborn-deep')

lower_idx = [e2-e1 for (e1, e2) in zip([fwi_paradise_confint[0],isi_paradise_confint[0],bui_paradise_confint[0],hdw_paradise_confint[0],fos_paradise_confint[0]],[FWI_paradise,ISI_paradise,BUI_paradise,HDW_paradise,FOS_paradise])]
upper_idx = [e1-e2 for (e1, e2) in zip([fwi_paradise_confint[1],isi_paradise_confint[1],bui_paradise_confint[1],hdw_paradise_confint[1],fos_paradise_confint[1]],[FWI_paradise,ISI_paradise,BUI_paradise,HDW_paradise,FOS_paradise])]

lower_met = [e2-e1 for (e1,e2) in zip([tmax_paradise_confint[0],vpd_paradise_confint[0],rh_paradise_confint[0],ws_paradise_confint[0]],[TMAX_paradise,VPD_paradise,RH_paradise,WS_paradise])]
upper_met = [e1-e2 for (e1,e2) in zip([tmax_paradise_confint[1],vpd_paradise_confint[1],rh_paradise_confint[1],ws_paradise_confint[1]],[TMAX_paradise,VPD_paradise,RH_paradise,WS_paradise])]

plt.figure(num=None, figsize=(17, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.9, wspace=0.3, hspace=0.3)
plt.rcParams.update({'font.size': 20})
ax1 = plt.subplot(1,1,1)
ax1.bar([6,7,8,9], [TMAX_paradise,VPD_paradise,RH_paradise,WS_paradise],yerr=[lower_met,upper_met], width=0.8,color='darkgrey')
ax1.bar([1,2,3,4,5], [FWI_paradise,ISI_paradise,BUI_paradise,HDW_paradise,FOS_paradise],yerr=[lower_idx,upper_idx],width=0.8,color='darkgrey')
ax1.set_ylim([-20,130])
ax1.set_ylabel('change in frequency')
ax2= ax1.twinx()
ax2.plot([0,10],[0,0],color='k')
#ax2.plot([0,10],[0.05,0.05],color='k',linestyle='--')
ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel('return period')
ax2.set_yticks([0.0,0.3333,1])
ax2.set_yticklabels(['1/20','1/15','1/10'])

plt.xticks([1,2,3,4,5,6,7,8,9],['FWI','ISI','BUI','HDW','FFWI','TA','VPD','RH','WS'])

plt.title('Paradise Fire',fontsize=24)
plt.box(False)

plotname = 'Paradise_barchart_PSAmean_deltaFrequency'
plt.savefig(plotname, dpi=300)

####### Bar chart Oregon ###########################
plt.style.use('seaborn-deep')

lower_idx = [e2-e1 for (e1, e2) in zip([fwi_oregon_confint[0],isi_oregon_confint[0],bui_oregon_confint[0],hdw_oregon_confint[0],fos_oregon_confint[0]],[FWI_oregon,ISI_oregon,BUI_oregon,HDW_oregon,FOS_oregon])]
upper_idx = [e1-e2 for (e1, e2) in zip([fwi_oregon_confint[1],isi_oregon_confint[1],bui_oregon_confint[1],hdw_oregon_confint[1],fos_oregon_confint[1]],[FWI_oregon,ISI_oregon,BUI_oregon,HDW_oregon,FOS_oregon])]

lower_met = [e2-e1 for (e1,e2) in zip([tmax_oregon_confint[0],vpd_oregon_confint[0],rh_oregon_confint[0],ws_oregon_confint[0]],[TMAX_oregon,VPD_oregon,RH_oregon,WS_oregon])]
upper_met = [e1-e2 for (e1,e2) in zip([tmax_oregon_confint[1],vpd_oregon_confint[1],rh_oregon_confint[1],ws_oregon_confint[1]],[TMAX_oregon,VPD_oregon,RH_oregon,WS_oregon])]

plt.figure(num=None, figsize=(17, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.9, wspace=0.3, hspace=0.3)
plt.rcParams.update({'font.size': 20})
ax1 = plt.subplot(1,1,1)
ax1.bar([6,7,8,9], [TMAX_oregon,VPD_oregon,RH_oregon,WS_oregon],yerr=[lower_met,upper_met], width=0.8,color='darkgrey')
ax1.bar([1,2,3,4,5], [FWI_oregon,ISI_oregon,BUI_oregon,HDW_oregon,FOS_oregon],yerr=[lower_idx,upper_idx],width=0.8,color='darkgrey')
ax1.set_ylim([-20,130])
ax1.set_ylabel('change in frequency')
ax2= ax1.twinx()
ax2.plot([0,10],[0,0],color='k')
#ax2.plot([0,10],[0.05,0.05],color='k',linestyle='--')
ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel('return period')
ax2.set_yticks([0.0,0.3333,1])
ax2.set_yticklabels(['1/20','1/15','1/10'])

plt.xticks([1,2,3,4,5,6,7,8,9],['FWI','ISI','BUI','HDW','FFWI','TA','VPD','RH','WS'])

plt.title('Oregon Fire',fontsize=24)
plt.box(False)

plotname = 'Oregon_barchart_PSAmean_deltaFrequency'
plt.savefig(plotname, dpi=300)

####### Bar chart southern CA ###########################
plt.style.use('seaborn-deep')

lower_idx = [e2-e1 for (e1, e2) in zip([fwi_malibu_confint[0],isi_malibu_confint[0],bui_malibu_confint[0],hdw_malibu_confint[0],fos_malibu_confint[0]],[FWI_malibu,ISI_malibu,BUI_malibu,HDW_malibu,FOS_malibu])]
upper_idx = [e1-e2 for (e1, e2) in zip([fwi_malibu_confint[1],isi_malibu_confint[1],bui_malibu_confint[1],hdw_malibu_confint[1],fos_malibu_confint[1]],[FWI_malibu,ISI_malibu,BUI_malibu,HDW_malibu,FOS_malibu])]

lower_met = [e2-e1 for (e1,e2) in zip([tmax_malibu_confint[0],vpd_malibu_confint[0],rh_malibu_confint[0],ws_malibu_confint[0]],[TMAX_malibu,VPD_malibu,RH_malibu,WS_malibu])]
upper_met = [e1-e2 for (e1,e2) in zip([tmax_malibu_confint[1],vpd_malibu_confint[1],rh_malibu_confint[1],ws_malibu_confint[1]],[TMAX_malibu,VPD_malibu,RH_malibu,WS_malibu])]

plt.figure(num=None, figsize=(17, 6), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.15, bottom=.1, right=.9, top=.9, wspace=0.3, hspace=0.3)
plt.rcParams.update({'font.size': 20})
ax1 = plt.subplot(1,1,1)
ax1.bar([6,7,8,9], [TMAX_malibu,VPD_malibu,RH_malibu,WS_malibu],yerr=[lower_met,upper_met], width=0.8,color='darkgrey')
ax1.bar([1,2,3,4,5], [FWI_malibu,ISI_malibu,BUI_malibu,HDW_malibu,FOS_malibu],yerr=[lower_idx,upper_idx],width=0.8,color='darkgrey')
ax1.set_ylim([-20,130])
ax1.set_ylabel('change in frequency')
ax2= ax1.twinx()
ax2.plot([0,10],[0,0],color='k')
#ax2.plot([0,10],[0.05,0.05],color='k',linestyle='--')
ax2.set_ylim([-0.2,1.3])
ax2.set_ylabel('return period')
ax2.set_yticks([0.0,0.3333,1])
ax2.set_yticklabels(['1/20','1/15','1/10'])

plt.xticks([1,2,3,4,5,6,7,8,9],['FWI','ISI','BUI','HDW','FFWI','TA','VPD','RH','WS'])

plt.title('Malibu Fire',fontsize=24)
plt.box(False)

plotname = 'Malibu_barchart_PSAmean_deltaFrequency'
plt.savefig(plotname, dpi=300)

plt.show()

