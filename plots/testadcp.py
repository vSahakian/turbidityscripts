#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 09:23:39 2021

@author: vjs
"""


## Test read in ADCP data that was downloaded



from netCDF4 import Dataset as netcdf_dataset

import numpy as np
import matplotlib.pyplot as plt
# from cartopy import config
# import cartopy.crs as ccrs
# import cartopy.feature as feature
# import cmocean


## Paths
adcpPath = '/Users/vjs/turbidites/observational/data/testADCP/netcdf/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20181130T000002Z_20181130T235952Z-binMapNone.nc'

subplotsize = (12,3)

##############
## Get data from the netcdf file
dataset = netcdf_dataset(adcpPath)

## Get the variables from the netcdf dataset object:

# Get all values of depth:
depth = dataset.variables['depth'][:]    

# Get all values of time:
time = dataset.variables['time'][:]

# Get velocity for beam 1, which is a 2D array:
beam1vel = dataset.variables['velocity_beam1'][:,:]
    
## Make time and depth 2d variables to plot data against:
TIME,DEPTH = np.meshgrid(time,depth)

## Plot:
plt.pcolormesh(TIME,DEPTH,beam1vel.T,shading='nearest',cmap='coolwarm',vmin=-0.02,vmax=0.02)
plt.colorbar()
plt.xlabel('Time')
plt.ylabel('Depth (m)')
plt.show()