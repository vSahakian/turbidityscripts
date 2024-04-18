#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 19:56:50 2022

@author: vjs
"""
## Plot scatter and doulb ehistogram of 2019 turbidity events at BACAX



import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import seaborn as sns
sns.set()
import obspy as obs



#############
## parameters and paths

adcp_netcdfpath = '/Users/vjs/turbidites/observational/data/BACAX_adcp_2019/netCDF/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20191224T000005Z_20191224T235955Z-binMapNone.nc'
turbidity_path = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed/reformated_BACAX_ntu_2019_template_notfiltered_2_15_1.1254_1_0_0.5_20191224T035330_signal_only.mseed'


## Output figure path:
adcp_path_pdf = '/Users/vjs/turbidites/observational/figs/pdf/adcpDec2019.pdf'
adcp_path_png = '/Users/vjs/turbidites/observational/figs/png/adcpDec2019.png'

####### Plot parameters
## Figure size:
adcp_figsize = (16,12)
colormap = 'Spectral_r'

##### Other parameters:

## Start and end time to plot, as datetime objects:
plot_starttime = datetime(2019,12,24,2,50,0)
plot_endtime = datetime(2019,12,24,16,50,0)


# %%
### Read it in:

###### SOME NOTES:
## To access variable, do: objectname.variables.values()
## Object name is: objectname.variables.values().name, or .long_name
#### SHORTNAME FOR VARIABLE TYPES NEEDED:
## time: time
## u: East
## v: North
## w: Up
## meanBackscatter: beam averaged corrected backscatter
## binmap_depth: water depth of final velocity measurement bins
## depth: water depth of beam-averaged measurement bins, tilt-corrected
## range: range from transducer
    
## Other params needed for function:
seconds_in_day = 24*60*60 # 24 hours, 60 minutes, 60 seconds

## Read in ADCP data:
adcp = nc4.Dataset(adcp_netcdfpath)
## And turbidity:
turbidity = obs.read(turbidity_path)


# %% TO FUNCTIONIZE:
## Input:
## adcp netcdf
## turbidity stream
## starttme to plot as datetime object
## endttme to plot as datetime object

## Get start time of netcdf file:
adcp_stime = datetime.strptime(adcp.time_coverage_start,'%Y%m%dT%H%M%SZ')
#etime = UTCDateTime(adcp.time_coverage_end)



## Get values to plot:
    
## Get turbidity times and data to plot:
turb_times_full = []
for i_time in turbidity[0].times('UTCDateTime'):
    turb_times_full.append(i_time.datetime)    
## And get where they're between the start and end time specified:    
turb_times_full = np.array(turb_times_full)
turb_range_indices = np.where((turb_times_full > plot_starttime) & (turb_times_full < plot_endtime))[0]

## GEt time and data to plot:
turb_times = turb_times_full[turb_range_indices]
turb_data = turbidity[0].data[turb_range_indices]
    

## Get time in seconds from the start time:
seconds_from_start = adcp['time'][:].data*seconds_in_day - adcp['time'][:].data[0]*seconds_in_day
## Make an array of UTC datetime objects that are the start time:
adcp_time_full = []
for i in range(len(seconds_from_start)):
    adcp_time_full.append(adcp_stime + timedelta(seconds=seconds_from_start[i]))
## Convert to array:
adcp_time_full = np.array(adcp_time_full).astype(datetime)    


###### GET DATA WITHIN TIMEFRAME TO PLOT ######
## Get locations of time array in netcdf to plot:
time_range_indices = np.where((adcp_time_full > plot_starttime) & (adcp_time_full < plot_endtime))[0]

## Get tme range for ADCP data to plot:
adcp_time = adcp_time_full[time_range_indices]



## Get other values:
east = adcp['u'][:][time_range_indices]
north = adcp['v'][:][time_range_indices]
up = adcp['w'][:][time_range_indices]
bscat = adcp['meanBackscatter'][:][time_range_indices]
depth = adcp['depth'][:]



## Plot order:
plotorder = [east,north,up,bscat]
labelorder = ['East (m/s)','North (m/s)','Up (m/s)','Backscatter (dB)']

## Make a meshgrid:
TIME,DEPTH = np.meshgrid(adcp_time,depth.data)

## Set up axes to plot:
figure, axes = plt.subplots(nrows=4,ncols=1,figsize=adcp_figsize)
    
## For ach axes, plot:
for i_property_ind, i_property in enumerate(plotorder):
    i_axes = axes[i_property_ind]
    i_pmesh = i_axes.pcolormesh(TIME,DEPTH,i_property.data.T,cmap=colormap)
    i_axes.invert_yaxis()
    
    ## Add colorbar:
    i_cb = plt.colorbar(i_pmesh,ax=i_axes)
    ## Set colorbar label:
    i_cb.set_label(labelorder[i_property_ind],size=16)
    ## colorbar tick param font size: 
    i_cb.ax.tick_params(labelsize=14)
    
    ## Format x axis labels:
    i_axes.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d,%H:%M'))
    i_axes.xaxis.set_tick_params(labelsize=14,rotation=0)
    
    ## Set y-label:
    i_axes.yaxis.set_tick_params(labelsize=14)
    i_axes.set_ylabel('Water Depth (m)',fontsize=16)
    
    ## At the end, for the last one (backscatter), add turbidity...
    if i_property_ind == len(plotorder)-1:
        turbax = i_axes.twinx()
        turbax.plot(turb_times,turb_data,color='black',linewidth=1,label='Turbidity (NTU)')

        ## also add x label:
        i_axes.set_xlabel('UTC Date Time', fontsize=16)
        #i_axis.set_xlabel('Month/Day/Hour in 2019',fontsize=16)
        
## Save figure
figure.savefig(adcp_path_pdf,transparent=True)
figure.savefig(adcp_path_png,transparent=True)


# %%

