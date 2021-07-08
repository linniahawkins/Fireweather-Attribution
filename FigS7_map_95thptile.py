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
import scipy.io
from pyfunctions import map_ptile


# TMAX
ptile = 95.0
cmap = 'RdYlBu_r'
in_dir = './data/'
in_file = os.path.join(in_dir+'TMAX_allsets_b827_domain_byFWI.txt')
out_file = './plots/FWI_95thptile_natforcing.png'

map_ptile(in_file,out_file,ptile,'Tmax (C) 95th percentile',0,50,cmap)

# RH
ptile = 5.0
cmap = 'BrBG'
in_dir = './data/'
in_file = os.path.join(in_dir+'RH_allsets_b827_domain_byFWI.txt')
out_file = './plots/RH_5thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'RH (%) 5th percentile',0,20,cmap)

# VPD
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'VPD_allsets_b827_domain_byFWI.txt')
out_file = './plots/VPD_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'VPD (hPa) 95th percentile',0,90,cmap)

# WS
ptile = 95.0
cmap = 'BuPu'
in_dir = './data/'
in_file = os.path.join(in_dir+'WIND_allsets_b827_domain_byFWI.txt')
out_file = './plots/WS_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'WS (m/s) 95th percentile',0,60,cmap)

# FWI
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'FWI_allsets_b827_domain_byFWI.txt')
out_file = './plots/FWI_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'FWI 95th percentile',0,220,cmap)

# ISI
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'ISI_allsets_b827_domain_byFWI.txt')
out_file = './plots/ISI_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'ISI 95th percentile',0,200,cmap)

# BUI
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'BUI_allsets_b827_domain_byFWI.txt')
out_file = './plots/BUI_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'BUI 95th percentile',0,1000,cmap)

# HDW
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'HDW_allsets_b827_domain_byFWI.txt')
out_file = './plots/HDW_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'HDW 95th percentile',0,3000,cmap)

# Fosberg
ptile = 95.0
cmap = 'YlOrBr'
in_dir = './data/'
in_file = os.path.join(in_dir+'FOS_allsets_b827_domain_byFWI.txt')
out_file = './plots/FOS_95thptile_fullforcing.png'

map_ptile(in_file,out_file,ptile,'Fosberg 95th percentile',0,120,cmap)
