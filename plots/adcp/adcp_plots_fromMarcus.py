#!/usr/bin/env python
# coding: utf-8

# ## ADCP Plots
# 
# Makes ADCP plots. Filters out background noise and scales colorbar by 90th percentile of max.
# 
# Input:
# - netCDF files, different frequencies and locations
# - use files from directory marcus_data
# 
# 
# Output:
# - ADCP plots from thesis

# In[35]:


import numpy as np
import matplotlib.pyplot as plt
import re
import os
import pandas as pd
import xarray as xr
import openpyxl
import netCDF4
from mpl_toolkits import mplot3d
from obspy.core.utcdatetime import UTCDateTime
from netCDF4 import Dataset as netcdf_dataset
import matplotlib.dates as mdates
import seaborn as sns
sns.set()


# In[143]:


#Paths
adcppath = '/Users/mjcu11/Downloads/search25175950/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20191223T000000Z_20191226T000000Z-Ensemble60s_binMapNone.nc'
ADCPPATH2MHZ_axis = '/Users/vjs/turbidites/marcus_data/AXIS/AXIS_ADCP_2MHz_Oct14-Oct28_2018/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc'
# ADCPPATH2MHZ_mideast = '/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_2MHz_Oct14-28-2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc'
# ADCPPATH2MHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_2MHz_Oct14-28_2018/BarkleyCanyonUpperSlope_UpperSlopeSouth_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc'

# ADCPPATH55KHZ_axis = '/Users/mjcu11/Documents/Val lab work/data/AXIS/AXIS_ADCP_55kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler55kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone_PLAN0.nc'

# ADCPPATH75KHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_75kHZ_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonUpperSlope_AcousticDopplerCurrentProfiler75kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'

# ADCPPATH150KHZ_mideast = '/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_150kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'

# ADCPPATH600KHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_600KHZ_Oct14-28_2018/search24702374/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'


# In[139]:


#ADCP IN XARRAY
ADCP55MHZ = xr.open_dataset('/Users/mjcu11/Downloads/search25175950/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20191223T000000Z_20191226T000000Z-Ensemble60s_binMapNone.nc')
ADCP55MHZ.transpose()

ADCP2MHZ_axis = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/AXIS/AXIS_ADCP_2MHz_Oct14-Oct28_2018/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc')
ADCP2MHZ_axis.transpose()

ADCP2MHZ_mideast = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_2MHz_Oct14-28-2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc')
ADCP2MHZ_mideast.transpose()

ADCP2MHZ_upperslope = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_2MHz_Oct14-28_2018/BarkleyCanyonUpperSlope_UpperSlopeSouth_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc')
ADCP2MHZ_upperslope.transpose()

ADCP55KHZ_axis = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/AXIS/AXIS_ADCP_55kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler55kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone_PLAN0.nc') 
ADCP55KHZ_axis.transpose()

ADCP75KHZ_upperslope = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_75kHZ_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonUpperSlope_AcousticDopplerCurrentProfiler75kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc')
ADCP75KHZ_upperslope.transpose()

ADCP150KHZ_mideast = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_150kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc')
ADCP150KHZ_mideast.transpose()

ADCP600KHZ_upperslope = xr.open_dataset('/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_600KHZ_Oct14-28_2018/search24702374/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc')
ADCP600KHZ_upperslope.transpose()


# In[140]:


#SET VARIABLES
#adcp
time55MHZ = ADCP55MHZ['time'].data
bin_depths55MHZ = ADCP55MHZ['binmap_depth'].data
east55MHZ = ADCP55MHZ['velocity_beam1'].data.T
north55MHZ = ADCP55MHZ['velocity_beam2'].data.T
up55MHZ = ADCP55MHZ['velocity_beam3'].data.T
rnge55MHZ = ADCP55MHZ['range'].data
u55MHZ = ADCP55MHZ['u'].data.T
v55MHZ = ADCP55MHZ['v'].data.T
w55MHZ = ADCP55MHZ['w'].data.T


time2MHZ_axis = ADCP2MHZ_axis['time'].data
bin_depths2MHZ_axis = ADCP2MHZ_axis['binmap_depth'].data
east2MHZ_axis = ADCP2MHZ_axis['velocity_beam1'].data.T
north2MHZ_axis = ADCP2MHZ_axis['velocity_beam2'].data.T
up2MHZ_axis = ADCP2MHZ_axis['velocity_beam3'].data.T
rnge2MHZ_axis = ADCP2MHZ_axis['range'].data
u2MHZ_axis = ADCP2MHZ_axis['u'].data.T
v2MHZ_axis = ADCP2MHZ_axis['v'].data.T
w2MHZ_axis = ADCP2MHZ_axis['w'].data.T

time2MHZ_mideast = ADCP2MHZ_mideast['time'].data
bin_depths2MHZ_mideast = ADCP2MHZ_mideast['binmap_depth'].data
east2MHZ_mideast = ADCP2MHZ_mideast['velocity_beam1'].data.T
north2MHZ_mideast = ADCP2MHZ_mideast['velocity_beam2'].data.T
up2MHZ_mideast = ADCP2MHZ_mideast['velocity_beam3'].data.T
rnge2MHZ_mideast = ADCP2MHZ_mideast['range'].data
u2MHZ_mideast = ADCP2MHZ_mideast['u'].data.T
v2MHZ_mideast = ADCP2MHZ_mideast['v'].data.T
w2MHZ_mideast = ADCP2MHZ_mideast['w'].data.T

time2MHZ_upperslope = ADCP2MHZ_upperslope['time'].data
bin_depths2MHZ_upperslope = ADCP2MHZ_upperslope['binmap_depth'].data
east2MHZ_upperslope = ADCP2MHZ_upperslope['velocity_beam1'].data.T
north2MHZ_upperslope = ADCP2MHZ_upperslope['velocity_beam2'].data.T
up2MHZ_upperslope = ADCP2MHZ_upperslope['velocity_beam3'].data.T
rnge2MHZ_upperslope = ADCP2MHZ_upperslope['range'].data
u2MHZ_upperslope = ADCP2MHZ_upperslope['u'].data.T
v2MHZ_upperslope = ADCP2MHZ_upperslope['v'].data.T
w2MHZ_upperslope = ADCP2MHZ_upperslope['w'].data.T

time55KHZ_axis = ADCP55KHZ_axis['time'].data
bin_depths55KHZ_axis = ADCP55KHZ_axis['binmap_depth'].data
east55KHZ_axis = ADCP55KHZ_axis['velocity_beam1'].data.T
north55KHZ_axis = ADCP55KHZ_axis['velocity_beam2'].data.T
up55KHZ_axis = ADCP55KHZ_axis['velocity_beam3'].data.T
rnge55KHZ_axis = ADCP55KHZ_axis['range'].data
u55KHZ_axis = ADCP55KHZ_axis['u'].data.T
v55KHZ_axis = ADCP55KHZ_axis['v'].data.T
w55KHZ_axis = ADCP55KHZ_axis['w'].data.T

time75KHZ_upperslope = ADCP75KHZ_upperslope['time'].data
bin_depths75KHZ_upperslope = ADCP75KHZ_upperslope['binmap_depth'].data
east75KHZ_upperslope = ADCP75KHZ_upperslope['velocity_beam1'].data.T
north75KHZ_upperslope = ADCP75KHZ_upperslope['velocity_beam2'].data.T
up75KHZ_upperslope = ADCP75KHZ_upperslope['velocity_beam3'].data.T
rnge75KHZ_upperslope = ADCP75KHZ_upperslope['range'].data
u75KHZ_upperslope = ADCP75KHZ_upperslope['u'].data.T
v75KHZ_upperslope = ADCP75KHZ_upperslope['v'].data.T
w75KHZ_upperslope = ADCP75KHZ_upperslope['w'].data.T

time150KHZ_mideast = ADCP150KHZ_mideast['time'].data
bin_depths150KHZ_mideast = ADCP150KHZ_mideast['binmap_depth'].data
east150KHZ_mideast = ADCP150KHZ_mideast['velocity_beam1'].data.T
north150KHZ_mideast = ADCP150KHZ_mideast['velocity_beam2'].data.T
up150KHZ_mideast = ADCP150KHZ_mideast['velocity_beam3'].data.T
rnge150KHZ_mideast = ADCP150KHZ_mideast['range'].data
u150KHZ_mideast = ADCP150KHZ_mideast['u'].data.T
v150KHZ_mideast = ADCP150KHZ_mideast['v'].data.T
w150KHZ_mideast = ADCP150KHZ_mideast['w'].data.T

time600KHZ_upperslope = ADCP600KHZ_upperslope['time'].data
bin_depths600KHZ_upperslope = ADCP600KHZ_upperslope['binmap_depth'].data
east600KHZ_upperslope = ADCP600KHZ_upperslope['velocity_beam1'].data.T
north600KHZ_upperslope = ADCP600KHZ_upperslope['velocity_beam2'].data.T
up600KHZ_upperslope = ADCP600KHZ_upperslope['velocity_beam3'].data.T
rnge600KHZ_upperslope = ADCP600KHZ_upperslope['range'].data
u600KHZ_upperslope = ADCP600KHZ_upperslope['u'].data.T
v600KHZ_upperslope = ADCP600KHZ_upperslope['v'].data.T
w600KHZ_upperslope = ADCP600KHZ_upperslope['w'].data.T


#Screen out surface noise velocities and set colorbar limit to 90th percentile of the data
lim_east55MHZ = float("%2.2f" % np.nanpercentile(east2MHZ , 90))
lim_north55MHZ = float("%2.2f" % np.nanpercentile(north2MHZ , 90))
lim_up55MHZ = float("%2.2f" % np.nanpercentile(up2MHZ , 90))
u_v_w55MHZ = max([lim_east2MHZ , lim_north2MHZ , lim_up2MHZ ])


lim_east2MHZ_axis = float("%2.2f" % np.nanpercentile(east2MHZ_axis , 90))
lim_north2MHZ_axis = float("%2.2f" % np.nanpercentile(north2MHZ_axis , 90))
lim_up2MHZ_axis = float("%2.2f" % np.nanpercentile(up2MHZ_axis , 90))
u_v_w2MHZ_axis = max([lim_east2MHZ_axis , lim_north2MHZ_axis , lim_up2MHZ_axis ])

lim_east2MHZ_mideast = float("%2.2f" % np.nanpercentile(east2MHZ_mideast , 90))
lim_north2MHZ_mideast = float("%2.2f" % np.nanpercentile(north2MHZ_mideast , 90))
lim_up2MHZ_mideast = float("%2.2f" % np.nanpercentile(up2MHZ_mideast , 90))
u_v_w2MHZ_mideast = max([lim_east2MHZ_mideast , lim_north2MHZ_mideast , lim_up2MHZ_mideast ])

lim_east2MHZ_upperslope = float("%2.2f" % np.nanpercentile(east2MHZ_upperslope , 90))
lim_north2MHZ_upperslope = float("%2.2f" % np.nanpercentile(north2MHZ_upperslope , 90))
lim_up2MHZ_upperslope = float("%2.2f" % np.nanpercentile(up2MHZ_upperslope , 90))
u_v_w2MHZ_upperslope = max([lim_east2MHZ_upperslope , lim_north2MHZ_upperslope , lim_up2MHZ_upperslope ])

lim_east55KHZ_axis = float("%2.2f" % np.nanpercentile(east55KHZ_axis, 90))
lim_north55KHZ_axis = float("%2.2f" % np.nanpercentile(north55KHZ_axis, 90))
lim_up55KHZ_axis = float("%2.2f" % np.nanpercentile(up55KHZ_axis, 90))
u_v_w55KHZ_axis = max([lim_east55KHZ_axis, lim_north55KHZ_axis, lim_up55KHZ_axis])

lim_east75KHZ_upperslope = float("%2.2f" % np.nanpercentile(east75KHZ_upperslope, 90))
lim_north75KHZ_upperslope = float("%2.2f" % np.nanpercentile(north75KHZ_upperslope, 90))
lim_up75KHZ_upperslope = float("%2.2f" % np.nanpercentile(up75KHZ_upperslope, 90))
u_v_w75KHZ_upperslope = max([lim_east75KHZ_upperslope, lim_north75KHZ_upperslope, lim_up75KHZ_upperslope])

lim_east150KHZ_mideast = float("%2.2f" % np.nanpercentile(east150KHZ_mideast, 90))
lim_north150KHZ_mideast = float("%2.2f" % np.nanpercentile(north150KHZ_mideast, 90))
lim_up150KHZ_mideast = float("%2.2f" % np.nanpercentile(up150KHZ_mideast, 90))
u_v_w150KHZ_mideast = max([lim_east150KHZ_mideast, lim_north150KHZ_mideast, lim_up150KHZ_mideast])

lim_east600KHZ_upperslope = float("%2.2f" % np.nanpercentile(east600KHZ_upperslope, 90))
lim_north600KHZ_upperslope = float("%2.2f" % np.nanpercentile(north600KHZ_upperslope, 90))
lim_up600KHZ_upperslope = float("%2.2f" % np.nanpercentile(up600KHZ_upperslope, 90))
u_v_w600KHZ_upperslope = max([lim_east600KHZ_upperslope, lim_north600KHZ_upperslope, lim_up600KHZ_upperslope])


# In[150]:


fig, axes = plt.subplots(3, sharex=True,)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon Axis ADCP 2 MHz",fontsize=24)

p0 = axes[0].pcolormesh(time55MHZ, rnge55MHZ, u55MHZ, cmap='seismic',vmin=-u_v_w55MHZ,vmax=u_v_w55MHZ)
#axes[0].set_ylim(0,800)
#axes[0].set_xlim('2018-10-22T00:05:39.000Z','2018-10-24T23:09:30.000Z')
axes[0].set_title('Eastward Velocity (u)')
#axes[0].annotate('three seismic events', xy=('2018-10-22-06',1),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',1.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
#axes[0].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#ax.annotate('ocean floor', xy=(6.28, 1), xytext=(10, 4),arrowprops=dict(facecolor='black', shrink=0.05,bbox=dict(boxstyle="round", alpha=0.1)))

#axes[0].invert_yaxis()
#axes[0].axvline('2019-12-23T20:56:23.302Z', color='yellow',linestyle='--',linewidth=3)
axes[0].axvline('2019-12-24T03:45:23.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time55MHZ, rnge55MHZ, v55MHZ, cmap='seismic',vmin=-u_v_w55MHZ,vmax=u_v_w55MHZ)
axes[1].set_title('Northward Velocity (v)')
#axes[1].set_ylim(0,800)
#axes[1].invert_yaxis()
axes[1].set_ylabel('Depth (m)')
#axes[1].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
axes[1].axvline('2019-12-24T03:45:23.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time55MHZ, rnge55MHZ, w55MHZ, cmap='seismic')
axes[2].set_title('Upward Velocity (w)')
#axes[2].set_ylim(0.1,0.7)
#axes[2].set_xlim('2018-10-18','2018-10-19')
#axes[2].invert_yaxis()
axes[2].axvline('2019-12-24T03:45:23.302Z', color='yellow',linestyle='--',linewidth=3)
axes[2].axvline('2019-12-25T03:50:23.302Z', color='yellow',linestyle='--',linewidth=3)
#axes[2].annotate('ocean floor', xy=('2018-10-27',0.12),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].set_xlabel('Time (Y-M-D)')

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[157]:


fig, axes = plt.subplots(nrows=3,sharex=True)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon MidEast ADCP 2 MHz",fontsize=24)

p0 = axes[0].pcolormesh(time2MHZ_mideast, rnge2MHZ_mideast, u2MHZ_mideast, cmap='seismic',vmin=-u_v_w2MHZ_mideast,vmax=u_v_w2MHZ_mideast)
#fig.colorbar(im, ax=ax0)
#ax0.set_ylim(0,800)
axes[0].set_title('Eastward Velocity (u)')
axes[0].annotate('three seismic events', xy=('2018-10-22-06',1),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',1.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
axes[0].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#axes[0].invert_yaxis()
axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time2MHZ_mideast, rnge2MHZ_mideast, v2MHZ_mideast, cmap='seismic',vmin=-u_v_w2MHZ_mideast,vmax=u_v_w2MHZ_mideast)
#fig.colorbar(im1, ax=ax1)
#ax1.set_ylim(0,800)
axes[1].set_title('Northward Velocity (v)')
axes[1].set_ylabel('Depth (m)')
axes[1].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#axes[1].invert_yaxis()
axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time2MHZ_mideast, rnge2MHZ_mideast, w2MHZ_mideast, cmap='seismic')
#fig.colorbar(im2, ax=ax2)
#ax2.set_ylim(800,1200)
axes[2].set_title('Upward Velocity (w)')
axes[2].set_xlabel('Time (Y-M-D)')
axes[2].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#axes[2].invert_yaxis()
axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[46]:


fig, axes = plt.subplots(nrows=3,sharex=True)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon Upperslope ADCP 2 MHz",fontsize=24)

p0 = axes[0].pcolormesh(time2MHZ_upperslope, rnge2MHZ_upperslope, u2MHZ_upperslope, cmap='seismic',vmin=-u_v_w2MHZ_upperslope,vmax=u_v_w2MHZ_upperslope)
#fig.colorbar(im, ax=ax0)
#ax0.set_ylim(0,800)
axes[0].set_title('Eastward Velocity (u)')
axes[0].annotate('three seismic events', xy=('2018-10-22-06',1),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',1.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
#axes[0].invert_yaxis()
axes[0].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time2MHZ_upperslope, rnge2MHZ_upperslope, v2MHZ_upperslope, cmap='seismic',vmin=-u_v_w2MHZ_upperslope,vmax=u_v_w2MHZ_upperslope)
#fig.colorbar(im1, ax=ax1)
#ax1.set_ylim(0,800)
axes[1].set_title('Northward Velocity (v)')
axes[1].set_ylabel('Depth (m)')
#axes[1].invert_yaxis()
axes[1].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time2MHZ_upperslope, rnge2MHZ_upperslope, w2MHZ_upperslope, cmap='seismic')
#fig.colorbar(im2, ax=ax2)
#ax2.set_ylim(800,1200)
axes[2].set_title('Upward Velocity (w)')
axes[2].set_xlabel('Time (Y-M-D)')
#axes[2].invert_yaxis()
axes[2].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[178]:


fig, axes = plt.subplots(3, sharex=True,)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon Axis ADCP 55kHZ",fontsize=24)

p0 = axes[0].pcolormesh(time55KHZ_axis, rnge55KHZ_axis, u55KHZ_axis, cmap='seismic',vmin=-u_v_w55KHZ_axis,vmax=u_v_w55KHZ_axis)
axes[0].set_ylim(0,900)
#axes[0].set_xlim('2018-10-22T00:05:39.000Z','2018-10-24T23:09:30.000Z')
axes[0].set_title('Eastward Velocity (u)')
axes[0].annotate('three seismic events', xy=('2018-10-22-06',800),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',900),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
axes[0].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#ax.annotate('ocean floor', xy=(6.28, 1), xytext=(10, 4),arrowprops=dict(facecolor='black', shrink=0.05,bbox=dict(boxstyle="round", alpha=0.1)))

#axes[0].invert_yaxis()
axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time55KHZ_axis, rnge55KHZ_axis, v55KHZ_axis, cmap='seismic',vmin=-u_v_w55KHZ_axis,vmax=u_v_w55KHZ_axis)
axes[1].set_title('Northward Velocity (v)')
axes[1].set_ylim(0,900)
#axes[1].invert_yaxis()
axes[1].set_ylabel('Depth (m)')
axes[1].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time55KHZ_axis, rnge55KHZ_axis, w55KHZ_axis, cmap='seismic')
axes[2].set_title('Upward Velocity (w)')
axes[2].set_ylim(0,900)
#axes[2].set_xlim('2018-10-18','2018-10-19')
#axes[2].invert_yaxis()
axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)
axes[2].annotate('ocean floor', xy=('2018-10-27',0.12),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].set_xlabel('Time (Y-M-D)')

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[188]:


fig, axes = plt.subplots(3, sharex=True,)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon MidEast ADCP 150kHZ",fontsize=24)

p0 = axes[0].pcolormesh(time150KHZ_mideast, rnge150KHZ_mideast, u150KHZ_mideast, cmap='seismic',vmin=-u_v_w150KHZ_mideast,vmax=u_v_w150KHZ_mideast)
#axes[0].set_ylim(0,900)
#axes[0].set_xlim('2018-10-22T00:05:39.000Z','2018-10-24T23:09:30.000Z')
axes[0].set_title('Eastward Velocity (u)')
axes[0].annotate('three seismic events', xy=('2018-10-22-06',200),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',250),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
axes[0].annotate('ocean floor', xy=('2018-10-27',10),bbox=dict(boxstyle="round",facecolor='white'))
#ax.annotate('ocean floor', xy=(6.28, 1), xytext=(10, 4),arrowprops=dict(facecolor='black', shrink=0.05,bbox=dict(boxstyle="round", alpha=0.1)))

#axes[0].invert_yaxis()
axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time150KHZ_mideast, rnge150KHZ_mideast, v150KHZ_mideast, cmap='seismic',vmin=-u_v_w150KHZ_mideast,vmax=u_v_w150KHZ_mideast)
axes[1].set_title('Northward Velocity (v)')
#axes[1].set_ylim(0,900)
#axes[1].invert_yaxis()
axes[1].set_ylabel('Depth (m)')
axes[1].annotate('ocean floor', xy=('2018-10-27',10),bbox=dict(boxstyle="round",facecolor='white'))
axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time150KHZ_mideast, rnge150KHZ_mideast, w150KHZ_mideast, cmap='seismic')
axes[2].set_title('Upward Velocity (w)')
#axes[2].set_ylim(0,900)
#axes[2].set_xlim('2018-10-18','2018-10-19')
#axes[2].invert_yaxis()
axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)
axes[2].annotate('ocean floor', xy=('2018-10-27',10),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].set_xlabel('Time (Y-M-D)')

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[202]:


fig, axes = plt.subplots(nrows=3,sharex=True)
fig.set_size_inches(20, 9)
fig.suptitle("Canyon Upperslope ADCP 75 kHZ",fontsize=24)

p0 = axes[0].pcolormesh(time75KHZ_upperslope, rnge75KHZ_upperslope, u75KHZ_upperslope, cmap='seismic',vmin=-u_v_w75KHZ_upperslope,vmax=u_v_w75KHZ_upperslope)
#fig.colorbar(im, ax=ax0)
#ax0.set_ylim(0,800)
axes[0].set_title('Eastward Velocity (u)')
axes[0].annotate('three seismic events', xy=('2018-10-22-06',350),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',410),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
#axes[0].invert_yaxis()
axes[0].annotate('ocean floor', xy=('2018-10-27',24),bbox=dict(boxstyle="round",facecolor='white'))
axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time75KHZ_upperslope, rnge75KHZ_upperslope, v75KHZ_upperslope, cmap='seismic',vmin=-u_v_w75KHZ_upperslope,vmax=u_v_w75KHZ_upperslope)
#fig.colorbar(im1, ax=ax1)
#ax1.set_ylim(0,800)
axes[1].set_title('Northward Velocity (v)')
axes[1].set_ylabel('Depth (m)')
#axes[1].invert_yaxis()
axes[1].annotate('ocean floor', xy=('2018-10-27',24),bbox=dict(boxstyle="round",facecolor='white'))
axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time75KHZ_upperslope, rnge75KHZ_upperslope, w75KHZ_upperslope, cmap='seismic')
#fig.colorbar(im2, ax=ax2)
#ax2.set_ylim(800,1200)
axes[2].set_title('Upward Velocity (w)')
axes[2].set_xlabel('Time (Y-M-D)')
#axes[2].invert_yaxis()
axes[2].annotate('ocean floor', xy=('2018-10-27',24),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[ ]:




