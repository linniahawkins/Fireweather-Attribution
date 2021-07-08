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

############ pyfunctions for Hawkins et al., GRL ##########

def bootstrap_confint_spatial(em1,em2, ptile, direction="above", c=[5, 95], bsn=1e5):
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
        sample_store[:,s] = np.divide(ful_occ,nat_occ)-1

    # now for each confidence interval find the frequency above/below ptile
    spatial_average = np.nanmean(sample_store,0)
    low = np.percentile(spatial_average, c[0],0)
    high = np.percentile(spatial_average, c[1],0)
    conf_inter = [low,high]
    
    return conf_inter

###########################################################


def bootstrap_confint_cell(em1,em2, ptile, direction="above", c=[0.05, 0.95], bsn=1e3):
    # ptile = percentile 
    # direction = "above" find frequency of occurance above ptile (or below)
    # c = confidence intervals (percentiles) to calculate
    # bsn = boot strap number, number of times to resample the distribution
    nat_data = em1.flatten()
    ful_data = em2.flatten()
    # create the store
    sample_store = np.zeros((int(bsn), 1), 'f')
    # define sample size
    #ssize = int(nat_data.shape[0] - 0.1*nat_data.shape[0])
    ssize = nat_data.shape[0]
    # do the resampling
    for s in range(0, int(bsn)):
        x1 = np.random.uniform(0,nat_data.shape[0],ssize).astype(int)
        x2 = np.random.uniform(0,nat_data.shape[0],ssize).astype(int)
        n_data = nat_data[x1]
        f_data = ful_data[x2]
        val = np.percentile(n_data,ptile)
        if direction == "above":
            nat_occ = float(len(n_data[n_data>val]))
            ful_occ = float(len(f_data[f_data>val]))
        elif direction == "below":
            nat_occ = float(len(n_data[n_data<val]))
            ful_occ = float(len(f_data[f_data<val]))
        else:
            print 'check direction input'
        sample_store[s] = np.divide(ful_occ,nat_occ)-1

    # now for each confidence interval find the frequency above/below ptile
    low = np.percentile(sample_store, c[0]*100)
    high = np.percentile(sample_store, c[1]*100)
    conf_inter = [low,high]
    return conf_inter


###########################################################

def bootstrap_domain(in_file1,in_file2,out_file,ptile,direction):
	# in_file1 : natural forcing data
	# in_file2 : full forcing data
	# out_file : save lowerbound of confint as
	# ptile    : percentile to bootstrap 
	# direction : above or below ptile 

	# get ocean mask
	oceanmask = np.genfromtxt('oceanmask.csv',delimiter=',')

    # readin data
	data = pd.read_csv(in_file1,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	natdata = data.values
	natptile = np.percentile(natdata,ptile,axis=1)
	natcount = [float(np.count_nonzero(natdata[cell,] > natptile[cell])) for cell in range(natdata.shape[0])]

	data = pd.read_csv(in_file2,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	fuldata = data.values
	fulcount = [float(np.count_nonzero(fuldata[cell,] > natptile[cell])) for cell in range(natdata.shape[0])]

	delta = np.divide(fulcount,natcount)
	fwi_domain = np.reshape(delta, (116,110))
	fwi_domain = fwi_domain+oceanmask

	######### find confidence interval by cell #########

	nat_domain = np.reshape(natdata, (116,110,2000))
	ful_domain = np.reshape(fuldata, (116,110,2000))

	fwi_confint = np.zeros([np.size(ful_domain,0),np.size(ful_domain,1),2])

	for lat in range(np.size(ful_domain,0)):
	    for lon in range(np.size(ful_domain,1)):
	        confint = bootstrap_confint_cell(nat_domain[lat,lon,:],ful_domain[lat,lon,:],ptile,direction=direction)
	        fwi_confint[lat,lon,0] = confint[0]
	        fwi_confint[lat,lon,1] = confint[1]

	np.savetxt(out_file,fwi_confint[:,:,0],delimiter=',',fmt='%1.6f')


###########################################################

def map_ratio(in_file1,in_file2,in_file3,out_file,ptile):
	# in_file1 : natural forcing simulations
	# in_file2 : full forcing simulations
	# in_file3 : lower bound of bootstrapped confidence interval (from bootstrap_domain)
	# out_file : save figure as
	# ptile    : percentile 

	# get ocean mask
	oceanmask = np.genfromtxt('oceanmask.csv',delimiter=',')

    # readin data
	data = pd.read_csv(in_file1,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	natdata = data.values
	natptile = np.percentile(natdata,ptile,axis=1)
	natcount = [float(np.count_nonzero(natdata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

	data = pd.read_csv(in_file2,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	fuldata = data.values
	fulcount = [float(np.count_nonzero(fuldata[cell,:] > natptile[cell])) for cell in range(natdata.shape[0])]

	delta = 100*(np.divide(fulcount,natcount)-1)
	fwi_domain = np.reshape(delta, (116,110))
	fwi_domain = fwi_domain+oceanmask

	######### readin confidence interval lower bound
	data = pd.read_csv(in_file3,index_col=None,header=None, sep=',',parse_dates=True, squeeze=True)                          
	lower = data.values
	mask = np.ma.masked_greater(lower, 0.01)
	mask = mask+oceanmask

	#### get lat lon
	dat = pd.read_csv('latlon.txt',index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	data = dat.transpose()
	lat = data.iloc[:,1].values
	lon = data.iloc[:,2].values
	lat2d = np.reshape(lat, (116, 110))
	lon2d = np.reshape(lon, (116, 110))

	# offset lat/lon
	file = 'latitude_longitude.nc'
	all_data=Dataset(file,'r')
	lat2d_winds=all_data.variables['global_latitude0']
	lon2d_winds=all_data.variables['global_longitude0']

	####### map ###################
	var = fwi_domain
	fig = plt.figure(num=None, figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
	plt.subplot(1,1,1)
	ax0 = plt.axes(projection=ccrs.Mercator())
	m = ax0.pcolormesh(lon2d,lat2d,var,transform=ccrs.PlateCarree(),cmap='RdGy_r',vmin=-200,vmax=200)

	b = ax0.pcolor(lon2d,lat2d, mask,transform=ccrs.PlateCarree(),hatch='xx', alpha=0.)

	ax0.set_extent([-125,-111,31.4,49.2])

	states = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '10m',
	                                edgecolor='black',
	                                facecolor='none')
	coast = cfeature.NaturalEarthFeature(category='physical', scale='10m', edgecolor='black',
	                        facecolor='none', name='coastline')
	ax0.add_feature(states)
	ax0.add_feature(coast)
	ax0.add_feature(cfeature.BORDERS, edgecolor='black')

	cbar = fig.colorbar(m)
	cbar.ax.tick_params(labelsize=20)
	cbar.ax.set_ylabel('change in frequency', fontsize=20)

	plt.savefig(out_file)


def map_ptile(in_file,out_file,ptile,varname,vmin,vmax,cmap):
	# in_file  : data
	# outfile  : save figure as
	# ptile    : percentile to map
	# varname  : variable name (i.e. FWI)
	# vmin     : lower bound color bar
	# vmax     : upper bound color bar
	# cmap     : color map

	#### get ocean mask
	oceanmask = np.genfromtxt('oceanmask.csv',delimiter=',')

	df = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	data = df.values
	datptile = np.percentile(data,ptile,axis=1)
	data_ptile = np.reshape(datptile, (116,110))+oceanmask

	#### get lat lon
	in_dir = './data/'
	in_file = os.path.join(in_dir+'latlon.txt')
	df = pd.read_csv(in_file,index_col=None,header=None, sep=' ',parse_dates=True, squeeze=True)                          
	data = df.transpose()
	lat = data.iloc[:,1].values
	lon = data.iloc[:,2].values
	lat2d = np.reshape(lat, (116, 110))
	lon2d = np.reshape(lon, (116, 110))

	#### Map 
	fig = plt.figure(num=None, figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
	plt.subplot(1,1,1)
	ax0 = plt.axes(projection=ccrs.Mercator())
	m = ax0.pcolormesh(lon2d,lat2d,data_ptile,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
	#draw continents and states
	ax0.set_extent([-125,-111,31.4,49.2])
	states = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines', '10m',
	                                edgecolor='black',
	                                facecolor='none')
	coast = cfeature.NaturalEarthFeature(category='physical', scale='10m', edgecolor='black',
	                        facecolor='none', name='coastline')
	ax0.add_feature(states)
	ax0.add_feature(coast)
	ax0.add_feature(cfeature.BORDERS, edgecolor='black')

	cbar = fig.colorbar(m)
	cbar.ax.tick_params(labelsize=20)
	cbar.ax.set_ylabel(varname, fontsize=20)

	plt.savefig(out_file)
