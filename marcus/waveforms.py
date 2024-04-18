#!/usr/bin/env python
# coding: utf-8

# ## Waveform plots
# 
# Downloads sample waveforms and makes time series plots and some preliminary maps

# In[ ]:


## http://ds.iris.edu/mda/NV/


# In[149]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
from glob import glob
import obspy as obs
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn.client import Client

# station
net = 'NV'
sta = 'BACME'

# time
stime = UTCDateTime("2018-10-22T053000.000")
etime = UTCDateTime("2018-10-22T063000.000")

# plot parameters:
figsize = (6,8)
figsize_long = (6,6)
  # limits
xlim = [stime.datetime,etime.datetime]
#bacme_ylim = [-138,152]
fraser_ylim = [-10.3,153]
bcheadturb_ylim = [0,2.25]
psalin_ylim = [34,34.06]
jma_ylim = [0,0.0042]
bcuturb_ylim = [0,25]
clay_ylim = [1033.285,1033.345]

## event times:
ev1 = UTCDateTime("2018-10-22T053939").datetime
ev2 = UTCDateTime("2018-10-22T061626").datetime
ev3 = UTCDateTime("2018-10-22T062248").datetime


# In[150]:


### First, download data for seismic station:
client = Client('IRIS')
bacme_st = client.get_waveforms(network=net,station=sta,location='*',channel='HNN,HNE,HNZ',attach_response=True,starttime=stime,endtime=etime)


# In[151]:


# remove response - get sensitivity for accelerometer since it's a flat gain:
bacme_sensitivity0 = bacme_st[0].stats.response.instrument_sensitivity.value
bacme_st[0].data = bacme_st[0].data / bacme_sensitivity0

## detrend it:
#bacme_st.detrend(type='linear')

## Get times for plotting, in utc format:
bacme_utctimes0 = bacme_st[0].times('utcdatetime')
bacme_datetimes0 = [bacme_utctimes0[i].datetime for i in range(len(bacme_st[0].times()))]
print(bacme_sensitivity0)


# In[152]:


# remove response - get sensitivity for accelerometer since it's a flat gain:
bacme_sensitivity1 = bacme_st[1].stats.response.instrument_sensitivity.value
bacme_st[1].data = bacme_st[1].data / bacme_sensitivity1

## detrend it:
#bacme_st.detrend(type='linear')

## Get times for plotting, in utc format:
bacme_utctimes1 = bacme_st[1].times('utcdatetime')
bacme_datetimes1 = [bacme_utctimes1[i].datetime for i in range(len(bacme_st[1].times()))]


# In[ ]:


# remove response - get sensitivity for accelerometer since it's a flat gain:
bacme_sensitivity2 = bacme_st[2].stats.response.instrument_sensitivity.value
bacme_st[2].data = bacme_st[2].data / bacme_sensitivity2

## detrend it:
#bacme_st.detrend(type='linear')

## Get times for plotting, in utc format:
bacme_utctimes2 = bacme_st[2].times('utcdatetime')
bacme_datetimes2 = [bacme_utctimes2[i].datetime for i in range(len(bacme_st[2].times()))]


# In[ ]:


baseline_shift = np.mean(bacme_st[0].data[0:10])
bacme_st[0].data = bacme_st[0].data - baseline_shift

baseline_shift = np.mean(bacme_st[1].data[0:10])
bacme_st[1].data = bacme_st[1].data - baseline_shift

baseline_shift = np.mean(bacme_st[2].data[0:10])
bacme_st[2].data = bacme_st[2].data - baseline_shift


# In[ ]:


pga0 = (abs(max(bacme_st[0].data))*0.10197)
print('pga east =',pga0)
pga1 = (abs(max(bacme_st[1].data))*0.10197)
print('pga north =',pga1)
pga2 = (abs(max(bacme_st[2].data))*0.10197)
print('pga up =',pga2)
pga=((pga0+pga1+pga2)/3)
print('pga =',pga)


# In[148]:


import matplotlib.pyplot as plt
import matplotlib.dates as md
####### Plot them:
## overall plot...
fig,axes = plt.subplots(3,1,sharex=True,constrained_layout=True,figsize=(10,8))
fig.suptitle("       CBC27",fontsize=34)
xformatter = md.DateFormatter('%H:%M')
xlocator = md.MinuteLocator(interval = 60)

## NV BACME
#East
axes[0].plot_date(bacme_datetimes0,bacme_st[0].data,'-',color='k',xdate=True,linewidth=0.5,label='NV.ENWF.HNE')
axes[0].axvline(ev1,color='red',linestyle='-',linewidth=1)
#axes[0].axvline(ev2,color='red',linestyle='-',linewidth=1)
#axes[0].axvline(ev3,color='red',linestyle='-',linewidth=1)
#axes[0].annotate('M6.5', xy=('5:40',-0.2),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('5:30',-0.2),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='blue'))
axes[0].grid('on')
axes[0].set_title('East',fontsize=27)
#axes[0].xaxis.set_major_locator(xlocator)
#axes[0].set_ylabel('Accel. (m/s/s)',labelpad=0,fontsize=27)
#axes[0].set_xlim(xlim)
#axes[0].set_ylim(bacme_ylim)
#axes[0].set_ylim(-0.03,0.03)
axes[0].legend(loc=1,fontsize=10,fancybox=True,facecolor='#e8e8ed',edgecolor='#3f3f40',shadow=True)

axes[1].plot_date(bacme_datetimes1,bacme_st[1].data,'-',color='k',xdate=True,linewidth=0.5,label='NV.NC89.HNN')
axes[1].axvline(ev1,color='red',linestyle='-',linewidth=1)
#axes[1].axvline(ev2,color='red',linestyle='-',linewidth=1)
#axes[1].axvline(ev3,color='red',linestyle='-',linewidth=1)
#axes[1].annotate('M6.3', xy=('0335',-0.2),bbox=dict(boxstyle="round",facecolor='yellow'),xytext=('03:55',-0.2),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-0.2",edgecolor='blue'))
axes[1].grid('on')
axes[1].set_title('North',fontsize=27)
#axes[1].set_xticklabels([])
#axes[1].set_ylabel('Accel. (m/s/s)',labelpad=0,fontsize=27)
#axes[1].set_xlim(xlim)
#axes[1].set_ylim(-0.08,0.08)
#axes[2].set_ylim(bacme_ylim)
axes[1].legend(loc=1,fontsize=10,fancybox=True,facecolor='#e8e8ed',edgecolor='#3f3f40',shadow=True)

axes[2].plot_date(bacme_datetimes2,bacme_st[2].data,'-',color='k',xdate=True,linewidth=0.5,label='NV.ENWF.HNZ')
axes[2].axvline(ev1,color='red',linestyle='-',linewidth=1)
axes[2].axvline(ev2,color='red',linestyle='-',linewidth=1)
axes[2].axvline(ev3,color='red',linestyle='-',linewidth=1)
axes[2].grid('on')
axes[2].set_title('ENWF',fontsize=27)
axes[2].set_ylabel('Accel. (m/s/s)',labelpad=0,fontsize=14)
#axes[2].set_xlabel('December 23rd, 2019',fontsize=27)
#axes[2].set_ylabel('Accel. (m/s/s)',labelpad=0,fontsize=27)
#axes[2].set_xlim(xlim)
axes[2].set_ylim(-0.015,0.015)
axes[2].legend(loc=1,fontsize=10,fancybox=True,facecolor='#e8e8ed',edgecolor='#3f3f40',shadow=True)

plt.gcf().axes[0].xaxis.set_major_formatter(xformatter)
plt.show()


# In[ ]:


#bacme_st[1].data[10]


# In[179]:


dict = {'datetime': bacme_datetimes0, 'east': bacme_st[0].data, 'north': bacme_st[1].data, 'up': bacme_st[2].data}  
       
df = pd.DataFrame(dict) 
    
# saving the dataframe 
df.to_csv('seis1.csv') 


# In[455]:


import os
import tarfile
from six.moves import urllib


# In[458]:


Endeavour_West_Flank_pga = 0.00290132985525
Barkley_mideast_pga = 0.00122489036923
Clayoquot_slope_pga = 0.00272776329758
Clayoquot_slope_bullseye_pga = 0.00278586031652
Cascadia_Basin_East_pga = 0.00420053063558
Cascadia_Basin_West_pga = 0.00168904713155
Barkley_canyon_pga = 0.00101615317537
Bamfield_pga = 9.63736426838e-05
Ucluelet_pga = 0.00011073307414
Endeavour_main_pga = 0.000215235756153
Tofino_pga = 0.000153020375463


# In[106]:


ds = pd.read_csv('/Users/mjcu11/Desktop/pga.csv',header=0)


# In[107]:


pga= ds.pga
lon= ds.lon
lat= ds.lat
loc= ds.location


# In[108]:


import matplotlib.pyplot as plt
ds.plot(kind="scatter", x='lon', y='lat', s=ds['pga']*100000,label="Peak Ground Acceleration (g)", c=ds['pga'], cmap=plt.get_cmap("jet"),
    colorbar=True, alpha=0.4, figsize=(10,7),)
plt.legend()
plt.show()


# In[110]:



import numpy as np
import matplotlib.image as mpimg
cascadia_img=mpimg.imread('/Users/mjcu11/Desktop/Screen Shot 2022-02-18 at 10.38.42 PM.png')
ax = ds.plot(kind="scatter", x='lon', y='lat', s=75,alpha=0.5,label="Peak Ground Acceleration (g)", c=ds['pga'], cmap=plt.get_cmap("cool"),
    colorbar=False, figsize=(10,7)
            )
plt.imshow(cascadia_img, extent=[-131, -124, 46, 50])
plt.ylabel("Latitude", fontsize=14)
plt.xlabel("Longitude", fontsize=14)

pga = ds["pga"]
tick_values = np.linspace(pga.min(), pga.max(), 11)
cbar = plt.colorbar()
#cbar.ax.set_yticklabels(["%d"%(round(v)) for v in tick_values], fontsize=14)
cbar.set_label('Peak Ground Acceleration (g)', fontsize=16)

plt.legend(fontsize=16)
plt.tight_layout()
plt.show()


# In[ ]:




