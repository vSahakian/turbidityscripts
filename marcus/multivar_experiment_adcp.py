#!/usr/bin/env python
# coding: utf-8

# ## V0, Multivariable and Experimenting with ADCP plots
# 
# Makes multivariable and ADCP plots at a certain station (or various stations)
# 
# Input: 
# - Turbidity
# - CTD
# - Oxygen
# - ADCP
# 
# Output: 
# - Multivariable plots
# - ADCP plots

# In[71]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import seaborn as sns
import datetime
import re
import os
import pandas as pd
import xarray as xr
import openpyxl
import netCDF4
from mpl_toolkits import mplot3d
from obspy.core.utcdatetime import UTCDateTime
from netCDF4 import Dataset as netcdf_dataset
import seaborn as sns


# In[99]:


#Paths
CHLOROTURB_axispath = pd.read_csv('/Users/mjcu11/Downloads/search25213177/BarkleyCanyon_BarkleyCanyonAxis_FluorometerTurbidity_Turbidity_20191223T000000Z_20191226T000000Z-NaN_clean_avg1minute.csv', skiprows=52)
CTD_axis = pd.read_csv('/Users/mjcu11/Documents/Val lab work/data/AXIS/AXIS_CTD_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonAxis_ConductivityTemperatureDepth_20181014T000000Z_20181028T000000Z-NaN_clean_avg10minute.csv', skiprows=52)
turbiditylong_axis= pd.read_csv('/Users/mjcu11/Downloads/search25084628/BarkleyCanyon_BarkleyCanyonAxis_FluorometerTurbidity_Turbidity_20180101T000000Z_20191231T000000Z-clean_avg15minute.csv')
O2_axispath = pd.read_csv('/Users/mjcu11/Downloads/search25213180/BarkleyCanyon_BarkleyCanyonAxis_OxygenSensor_20191223T000000Z_20191226T000000Z-NaN_clean_avg1minute.csv', skiprows=52)
O2_mideastpath = pd.read_csv('/Users/mjcu11/Downloads/search25213243/BarkleyCanyon_BarkleyCanyonMid-East_OxygenSensor_20191223T000000Z_20191226T000000Z-NaN_clean_avg1minute.csv',skiprows=52)
O2_hydratespath = pd.read_csv('/Users/mjcu11/Downloads/search25213209/BarkleyCanyon_BarkleyCanyonHydrates_OxygenSensor_20191223T000000Z_20191226T000000Z-NaN_clean_avg1minute.csv',skiprows=52)


# In[103]:


#variables
turb_axis = np.array(CHLOROTURB_axispath[' "Turbidity (NTU)"'])
time_axis = np.array(CHLOROTURB_axispath['"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"'])
O2_axis = np.array(O2_axispath[' "Oxygen Concentration Corrected (ml/l)"'])
temp_axis = np.array(O2_axispath[' "Temperature (C)"'])

O2_mideast = np.array(O2_mideastpath[' "Oxygen Concentration Corrected (ml/l)"'])
temp_mideast = np.array(O2_mideastpath[' "Temperature (C)"'])
time_mideast = np.array(O2_mideastpath['#"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"'])

O2_hydrates = np.array(O2_hydratespath[' "Oxygen Concentration Corrected (ml/l)"'])
temp_hydrates = np.array(O2_hydratespath[' "Temperature (C)"'])
time_hydrates = np.array(O2_hydratespath['#"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"'])


# In[72]:


sns.set()


# In[91]:


#chloroturb
fig, axes = plt.subplots(3,sharex=True,constrained_layout=True)
fig.suptitle("Barkley Axis 12/23 -- 12/26/2019",fontsize=24)
fig.set_size_inches(20, 9)

p0 = axes[0].plot(time_axis, turb_axis)
#axes[0].set_xlim('2018-10-14T05:00:00.302Z','2018-10-28T23:00:00.302Z')
axes[0].set_ylabel('Turbidity (NTU)',fontsize=20)
#axes[0].set_xticks(np.arange(0, 2016, step=1000))
#axes[0].set_xticks('2018-10-22T05:25:00.000Z', '2018-10-22T06:55:00.000Z')
#axes[0].set_xlim(1206, 1349)
#axes[0].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1243,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1245,color='red',linestyle='-',linewidth=1)
#axes[0].set_title('2018-10-14 -- 2018-10-28',fontsize=18)

p1 = axes[1].plot(time, O2_axis)
#axes[1].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:16:26.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:22:48.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].set_xlim('2018-10-22T05:00:00.302Z','2018-10-28T23:00:00.302Z')
#axes[1].set_xlim('2018-10-14T00:05:00.000Z','2018-10-27T23:55:00.000Z')
#axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%M-%D'))
axes[1].set_ylabel('O2 Conc. (ml/l)',fontsize= 20)
#axes[1].set_xlabel('Datetime',fontsize=18)

p2 = axes[2].plot(time, temp_axis)
#axes[1].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:16:26.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:22:48.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].set_xlim('2018-10-22T05:00:00.302Z','2018-10-28T23:00:00.302Z')
#axes[1].set_xlim('2018-10-14T00:05:00.000Z','2018-10-27T23:55:00.000Z')
axes[2].xaxis.set_major_locator(mdates.DayLocator(interval=300))
#axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%M-%D'))
axes[2].set_ylabel('Temp. (C)',fontsize= 20)
#axes[1].set_xlabel('Datetime',fontsize=18)

plt.gcf().autofmt_xdate()
plt.show()


# In[104]:


fig, axes = plt.subplots(2,sharex=True,constrained_layout=True)
fig.suptitle("Barkley Mideast 12/23 -- 12/26/2019",fontsize=24)
fig.set_size_inches(20, 9)

p0 = axes[0].plot(time_mideast, O2_hydrates)
#axes[0].set_xlim('2018-10-14T05:00:00.302Z','2018-10-28T23:00:00.302Z')
axes[0].set_ylabel('O2 Conc. (ml/l)',fontsize=20)
#axes[0].set_xticks(np.arange(0, 2016, step=1000))
#axes[0].set_xticks('2018-10-22T05:25:00.000Z', '2018-10-22T06:55:00.000Z')
#axes[0].set_xlim(1206, 1349)
#axes[0].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1243,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1245,color='red',linestyle='-',linewidth=1)
#axes[0].set_title('2018-10-14 -- 2018-10-28',fontsize=18)

p1 = axes[1].plot(time_hydrates, temp_hydrates)
#axes[1].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:16:26.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:22:48.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].set_xlim('2018-10-22T05:00:00.302Z','2018-10-28T23:00:00.302Z')
#axes[1].set_xlim('2018-10-14T00:05:00.000Z','2018-10-27T23:55:00.000Z')
#axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%M-%D'))
axes[1].xaxis.set_major_locator(mdates.DayLocator(interval=300))
axes[1].set_ylabel('Temp (C)',fontsize= 20)
#axes[1].set_xlabel('Datetime',fontsize=18)

plt.gcf().autofmt_xdate()
plt.show()


# In[105]:


fig, axes = plt.subplots(2,sharex=True,constrained_layout=True)
fig.suptitle("Barkley Hydrates 12/23 -- 12/26/2019",fontsize=24)
fig.set_size_inches(20, 9)

p0 = axes[0].plot(time_hydrates, O2_mideast)
#axes[0].set_xlim('2018-10-14T05:00:00.302Z','2018-10-28T23:00:00.302Z')
axes[0].set_ylabel('O2 Conc. (ml/l)',fontsize=20)
#axes[0].set_xticks(np.arange(0, 2016, step=1000))
#axes[0].set_xticks('2018-10-22T05:25:00.000Z', '2018-10-22T06:55:00.000Z')
#axes[0].set_xlim(1206, 1349)
#axes[0].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1243,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=0.5)
#axes[0].axvline(1244,color='red',linestyle='-',linewidth=1)
#axes[0].axvline(1245,color='red',linestyle='-',linewidth=1)
#axes[0].set_title('2018-10-14 -- 2018-10-28',fontsize=18)

p1 = axes[1].plot(time_mideast, temp_mideast)
#axes[1].axvline('2019-12-23T20:56:23.000Z', color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:16:26.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].axvline('2018-10-22T06:22:48.000Z',color='red',linestyle='-',linewidth=1)
#axes[1].set_xlim('2018-10-22T05:00:00.302Z','2018-10-28T23:00:00.302Z')
#axes[1].set_xlim('2018-10-14T00:05:00.000Z','2018-10-27T23:55:00.000Z')
#axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%M-%D'))
axes[1].xaxis.set_major_locator(mdates.DayLocator(interval=300))
axes[1].set_ylabel('Temp (C)',fontsize= 20)
#axes[1].set_xlabel('Datetime',fontsize=18)

plt.gcf().autofmt_xdate()
plt.show()


# In[119]:


#Paths
ADCPPATH2MHZ_axis = '/Users/mjcu11/Downloads/search25213578/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler55kHz_20191224T000000Z_20191224T010012Z-binMapNearest_PLAN0.nc'
ADCPPATH2MHZ_mideast = '/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_2MHz_Oct14-28-2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc'
ADCPPATH2MHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_2MHz_Oct14-28_2018/BarkleyCanyonUpperSlope_UpperSlopeSouth_AcousticDopplerCurrentProfiler2MHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone.nc'

ADCPPATH55KHZ_axis = '/Users/mjcu11/Documents/Val lab work/data/AXIS/AXIS_ADCP_55kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler55kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNone_PLAN0.nc'

ADCPPATH75KHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_75kHZ_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonUpperSlope_AcousticDopplerCurrentProfiler75kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'

ADCPPATH150KHZ_mideast = '/Users/mjcu11/Documents/Val lab work/data/MIDEAST/MIDEAST_ADCP_150kHz_Oct14-28_2018/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'

ADCPPATH600KHZ_upperslope = '/Users/mjcu11/Documents/Val lab work/data/UPPERSLOPE/UPPERSLOPE_ADCP_600KHZ_Oct14-28_2018/search24702374/BarkleyCanyon_BarkleyCanyonMid-East_AcousticDopplerCurrentProfiler150kHz_20181014T000000Z_20181028T000000Z-Ensemble600s_binMapNearest_3beamOn.nc'


# In[191]:


#ADCP IN XARRAY
ADCP2MHZ_axis = xr.open_dataset('/Users/mjcu11/Downloads/search25213639/BarkleyCanyon_BarkleyCanyonAxis_AcousticDopplerCurrentProfiler55kHz_20191223T000000Z_20191226T000000Z-Ensemble600s_binMapNone_PLAN0.nc')
ADCP2MHZ_axis.transpose()


# In[192]:


#SET VARIABLES
#adcp
time2MHZ_axis = ADCP2MHZ_axis['time'].data
bin_depths2MHZ_axis = ADCP2MHZ_axis['binmap_depth'].data
east2MHZ_axis = ADCP2MHZ_axis['velocity_beam1'].data.T
north2MHZ_axis = ADCP2MHZ_axis['velocity_beam2'].data.T
up2MHZ_axis = ADCP2MHZ_axis['velocity_beam3'].data.T
rnge2MHZ_axis = ADCP2MHZ_axis['range'].data
u2MHZ_axis = ADCP2MHZ_axis['u'].data.T
v2MHZ_axis = ADCP2MHZ_axis['v'].data.T
w2MHZ_axis = ADCP2MHZ_axis['w'].data.T







#Screen out surface noise velocities and set colorbar limit to 90th percentile of the data
lim_east2MHZ_axis = float("%2.2f" % np.nanpercentile(east2MHZ_axis , 90))
lim_north2MHZ_axis = float("%2.2f" % np.nanpercentile(north2MHZ_axis , 90))
lim_up2MHZ_axis = float("%2.2f" % np.nanpercentile(up2MHZ_axis , 90))
u_v_w2MHZ_axis = max([lim_east2MHZ_axis , lim_north2MHZ_axis , lim_up2MHZ_axis ])



# In[193]:



#Plots
#adcp
fig, axes = plt.subplots(3, sharex=True,)
fig.set_size_inches(20, 9)
fig.suptitle("Barkley Mideast ADCP 150 kHz (12/24/2019)",fontsize=24)

p0 = axes[0].pcolormesh(time2MHZ_axis, rnge2MHZ_axis, u2MHZ_axis, cmap='seismic',vmin=-u_v_w2MHZ_axis,vmax=u_v_w2MHZ_axis)
#axes[0].set_ylim(0,800)
#axes[0].set_xlim('2018-10-22T00:05:39.000Z','2018-10-24T23:09:30.000Z')
axes[0].set_title('Eastward Velocity (u)')
#axes[0].annotate('three seismic events', xy=('2018-10-22-06',1),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('2018-10-23',1.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='yellow'))
#axes[0].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#ax.annotate('ocean floor', xy=(6.28, 1), xytext=(10, 4),arrowprops=dict(facecolor='black', shrink=0.05,bbox=dict(boxstyle="round", alpha=0.1)))

axes[0].invert_yaxis()
#axes[0].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p1 = axes[1].pcolormesh(time2MHZ_axis, rnge2MHZ_axis, v2MHZ_axis, cmap='seismic',vmin=-u_v_w2MHZ_axis,vmax=u_v_w2MHZ_axis)
axes[1].set_title('Northward Velocity (v)')
#axes[1].set_ylim(0,800)
axes[1].invert_yaxis()
axes[1].set_ylabel('Depth (m)')
#axes[1].annotate('ocean floor', xy=('2018-10-27',0.15),bbox=dict(boxstyle="round",facecolor='white'))
#axes[1].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)

p2 = axes[2].pcolormesh(time2MHZ_axis, rnge2MHZ_axis, w2MHZ_axis, cmap='seismic')
axes[2].set_title('Upward Velocity (w)')
#axes[2].set_ylim(0.1,0.7)
#axes[2].set_xlim('2018-10-18','2018-10-19')
axes[2].invert_yaxis()
#axes[2].axvline('2018-10-22T05:39:00.302Z', color='yellow',linestyle='--',linewidth=3)
#axes[2].annotate('ocean floor', xy=('2018-10-27',0.12),bbox=dict(boxstyle="round",facecolor='white'))
axes[2].set_xlabel('Time (Y-M-D)')

fig.colorbar(p0,ax=axes.ravel().tolist(), label='Water Velocity m/s')

plt.show()


# In[132]:


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


# In[158]:


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


# In[3]:


#open multiple netcdf files and combine them into one netcdf
CTD_axis = xr.open_mfdataset('/Users/mjcu11/Downloads/search25104970 2/BarkleyCanyonUpperSlope*.nc',combine = 'by_coords')
CTD_axis.to_netcdf('CTD_combined.nc')

#turb_axis = xr.open_mfdataset('/Users/mjcu11/Downloads/search25135606/BarkleyCanyonUpperSlope_UpperSlopeVerticalProfiler_FluorometerTurbidity_20181022T050000Z_20181022T065959Z-clean.nc')


# In[335]:


#reopen single netcdf and continue as convert to pandas dataframe
ds = xr.open_mfdataset('CTD_combined.nc')
df = ds.to_dataframe()
#ds = xr.open_mfdataset('/Users/mjcu11/Downloads/search25135606/BarkleyCanyonUpperSlope_UpperSlopeVerticalProfiler_FluorometerTurbidity_20181022T050000Z_20181022T065959Z-clean.nc')
#df = ds.to_dataframe()
df


# In[336]:


#create variables
time_axis = CTD_axis['time'].data
depth_axis = CTD_axis['depth'].data
chloro_axis = CTD_axis['chlorophyll'].data
turb_axis = CTD_axis['turbidityntu'].data


# In[337]:


#add column
df['month'] = 'OCT'
df['year'] = '2018'
df['site ID'] = 'upperslope'
df['cruise'] = 'barkley'
df['lat'] = '48.427713'
df['lon'] = '-126.17408'
df['date'] = '10/22/18'
df['station'] = 1
df['time'] = time_axis


# In[338]:


# to remove column 
#df.drop("chlorophyll_qaqcFlags", axis=1, inplace=True)
#df.drop("depth_qaqcFlags", axis=1, inplace=True)
#df.drop("turbidityntu_qaqcFlags", axis=1, inplace=True)
#df.drop("pitch_qaqcFlags", axis=1, inplace=True)
#df.drop("roll_qaqcFlags", axis=1, inplace=True)

#to change name of column
df = df.rename(columns={"time":"time_ISO8601"})
df
#df.to_csv('turbshort.csv',date_format='%Y-%m-%dT%H:%M:%S.%f%z')


# In[339]:


df.to_csv('ctdnetcdf.csv',date_format='%Y-%m-%dT%H:%M:%S.%f%z')
#df.to_hdf('ctdnetcdf.hdf', key='stage', mode='w')


# In[321]:


ds= ds.from_dataframe(df)
ds


# In[275]:


ds= ds.rename({'time': 'time_ISO8601'})


# In[276]:


ds


# In[322]:


#export dataframe to netcdf file
import netCDF4 as nc
import numpy as np

turbnc =ds.to_netcdf(path='ctdnetcdf.nc', mode='w', format='NETCDF4', group=None, engine=None, encoding=None, unlimited_dims=None, compute=True, invalid_netcdf=False,)


# In[ ]:




