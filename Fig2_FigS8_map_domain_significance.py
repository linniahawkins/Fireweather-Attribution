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
from pyfunctions import map_ratio, bootstrap_domain

########## bootstrap domain ###################
 
# FWI
ptile = 95.0
in_dir = '../data/'
in_file1 = os.path.join(in_dir+'FWI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'fwi_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

# ISI
ptile = 95.0
in_file1 = os.path.join(in_dir+'ISI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'isi_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

# BUI
ptile = 95.0
in_file1 = os.path.join(in_dir+'BUI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'bui_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

# HDW
ptile = 95.0
in_file1 = os.path.join(in_dir+'HDW_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'hdw_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

# Fosberg
ptile = 95.0
in_file1 = os.path.join(in_dir+'FOS_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'fos_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")


# TMAX
ptile = 95.0
in_file1 = os.path.join(in_dir+'TMAX_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'tmax_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")


# RH
ptile = 5.0
in_file1 = os.path.join(in_dir+'RH_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'rh_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"below")

# VPD
ptile = 95.0
in_file1 = os.path.join(in_dir+'VPD_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'vpd_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

# WIND SPEED
ptile = 95.0
in_file1 = os.path.join(in_dir+'WIND_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
out_file = os.path.join(in_dir+'wind_domain_confint_lower.txt')

bootstrap_domain(in_file1,in_file2,out_file,ptile,"above")

#####################################################################
###### map ratio of frequency of Full>Natural 95th to Natural>Natural 95th
#####################################################################


############ map fwi  ##########
ptile = 95.0
in_dir = '../data/'
in_file1 = os.path.join(in_dir+'FWI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'fwi_domain_confint_lower.txt')
out_file = "fwi_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map isi  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'ISI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'isi_domain_confint_lower.txt')
out_file = "isi_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map bui  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'BUI_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'bui_domain_confint_lower.txt')
out_file = "bui_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map hdw  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'HDW_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'hdw_domain_confint_lower.txt')
out_file = "hdw_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map fosberg  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'FOS_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'fos_domain_confint_lower.txt')
out_file = "fos_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)


############ map tmax  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'TMAX_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'tmax_domain_confint_lower.txt')
out_file = "tmax_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)


############ map rh  ##########
ptile = 5.0
in_file1 = os.path.join(in_dir+'RH_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'rh_domain_confint_lower.txt')
out_file = "rh_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map vpd  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'VPD_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'vpd_domain_confint_lower.txt')
out_file = "vpd_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

############ map windspeed  ##########
ptile = 95.0
in_file1 = os.path.join(in_dir+'WIND_allsets_b840_domain_byFWI.txt')
in_file2 = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
in_file3 = os.path.join(in_dir+'wind_domain_confint_lower.txt')
out_file = "wind_domain_bootstrap_significance_4.png"

map_ratio(in_file1,in_file2,in_file3,out_file,ptile)

