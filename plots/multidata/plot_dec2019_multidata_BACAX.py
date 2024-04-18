#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 19:56:50 2022

@author: vjs
"""
## Plot scatter and doulb ehistogram of 2019 turbidity events at BACAX



import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.dates import DateFormatter
from datetime import datetime, timedelta
import numpy as np

import matplotlib.dates as mdates
from glob import glob


#############
## parameters and paths

## Paths to all files for approx. time range:
turbidity_paths = sorted(glob('/Users/vjs/turbidites/observational/data/barkleycanyon/BACAX/turbidityntu/BarkleyCanyon_BarkleyCanyonAxis_variables_TurbidityNTU_2019122*_2019122*-clean_avg1minute.csv'))
temp_paths = sorted(glob('/Users/vjs/turbidites/observational/data/barkleycanyon/BACAX/seawatertemperature/BarkleyCanyon_BarkleyCanyonAxis_variables_SeaWaterTemperature_2019122*_2019122*-clean_avg1minute.csv'))
chlorophyll_paths = sorted(glob('/Users/vjs/turbidites/observational/data/barkleycanyon/BACAX/chlorophyll/BarkleyCanyon_BarkleyCanyonAxis_variables_Chlorophyll_2019122*_2019122*-clean_avg1minute.csv'))
oxygen_paths =  sorted(glob('/Users/vjs/turbidites/observational/data/barkleycanyon/BACAX/oxygen/BarkleyCanyon_BarkleyCanyonAxis_variables_Oxygen_2019122*_2019122*-clean_avg1minute.csv'))

events2019path = '/Users/vjs/turbidites/observational/data/templates/events_turbidity_BACAX2019.csv'

## Output figure path:
main_multidataplot_path_pdf = '/Users/vjs/turbidites/observational/figs/pdf/multidataDec2019_main.pdf'
main_multidataplot_path_png = '/Users/vjs/turbidites/observational/figs/png/multidataDec2019_main.png'

## Plotting params
figure_dims = (16,13)
# xmin = datetime(2019,12,23)
# xmax = datetime(2019,12,24,23,59,00)
xmin = datetime(2019,12,23,16,0,0)
xmax = datetime(2019,12,24,8,00,00)

## Earthquake times:
eq_path = '/Users/vjs/turbidites/observational/data/earthquakes/events_dec23_25_2019.csv'
    
    
# %%
## Read in turbidity events identified by STA/LTA:
events2019all = pd.read_csv(events2019path,names=['Time','Amp','Review'],skiprows=1)
## Get datetime column:
events2019all['FormattedTime'] = pd.to_datetime(events2019all.Time.astype('datetime64[ms]'))
## And take only vetted ones:
events2019 = events2019all.loc[events2019all.Review == 1]

## Also read in the earthquakes:
earthquakes_all = pd.read_csv(eq_path)
earthquakes_all['FormattedTime'] = pd.to_datetime(earthquakes_all.time.astype('datetime64[ms]'))
## Get larger than M6 events
earthquakes = earthquakes_all.loc[earthquakes_all.mag >= 6]

#########

### Read it in and put each value into one dataframe:
turbidity_df = pd.DataFrame(columns=['datetime_string','turbidityntu','QC','turbcount'])
for i_turbiditypath in turbidity_paths:
    i_data = pd.read_csv(i_turbiditypath,names=['datetime_string','turbidityntu','QC','turbcount'],comment='#')
    ## append:
    turbidity_df = pd.concat([turbidity_df,i_data],ignore_index=True)
## Get datetime not as a string:
turbidity_df['datetime'] = pd.to_datetime(turbidity_df['datetime_string'].astype('datetime64[ms]'))

## For temperature:    
temperature_df = pd.DataFrame(columns=['datetime_string','seawatertemperature','QC','turbcount'])
for i_temperaturepath in temp_paths:
    i_data = pd.read_csv(i_temperaturepath,names=['datetime_string','seawatertemperature','QC','turbcount'],comment='#')
    ## append:
    temperature_df = pd.concat([temperature_df,i_data],ignore_index=True)
## Get datetime not as a string:
temperature_df['datetime'] = pd.to_datetime(temperature_df['datetime_string'].astype('datetime64[ms]'))
    
## For chlorophyll:
chlorophyll_df = pd.DataFrame(columns=['datetime_string','chlorophyll','QC','turbcount'])
for i_chlorophyllpath in chlorophyll_paths:
    i_data = pd.read_csv(i_chlorophyllpath,names=['datetime_string','chlorophyll','QC','turbcount'],comment='#')
    ## append:
    chlorophyll_df = pd.concat([chlorophyll_df,i_data],ignore_index=True)
## Get datetime not as a string:
chlorophyll_df['datetime'] = pd.to_datetime(chlorophyll_df['datetime_string'].astype('datetime64[ms]'))

## For oxygen:
oxygen_df = pd.DataFrame(columns=['datetime_string','oxygen','QC','turbcount'])
for i_oxygenpath in oxygen_paths:
    i_data = pd.read_csv(i_oxygenpath,names=['datetime_string','oxygen','QC','turbcount'],comment='#')
    ## append:
    oxygen_df = pd.concat([oxygen_df,i_data],ignore_index=True)
## Get datetime not as a string:
oxygen_df['datetime'] = pd.to_datetime(oxygen_df['datetime_string'].astype('datetime64[ms]'))


# %%
## Make multi sensor plot
multifig,multiaxes = plt.subplots(nrows=4,ncols=1,figsize=figure_dims,sharex=True)

dataframe_ordered = [oxygen_df,chlorophyll_df,temperature_df,turbidity_df]
column_ordered = ['oxygen','chlorophyll','seawatertemperature','turbidityntu']
ylabel_ordered = ['Oxygen (ml/L)','Chlorophyll (ug/L)','Sea Water \nTemperature (C)','Turbidity (NTU)']

## For each data type, plot:
for i_axis,i_dataframe,i_column,i_ylabel in zip(multiaxes,dataframe_ordered,column_ordered,ylabel_ordered):
    ## get axis
    # i_axis = multiaxes[i_axis]
    ## scatter plot:
    i_axis.plot(i_dataframe.datetime,i_dataframe[i_column],linewidth=1.6,color='#056e6c')
    
    ## set label
    i_axis.set_ylabel(i_ylabel,fontsize=16)
    
    ## Add lines for events identified:
    i_datamin = i_dataframe.loc[((i_dataframe.datetime>xmin) & (i_dataframe.datetime<xmax))][i_column].min()
    i_datamax = i_dataframe.loc[((i_dataframe.datetime>xmin) & (i_dataframe.datetime<xmax))][i_column].max()
    i_datadiff = i_datamax - i_datamin
    
    i_ymin = i_datamin - i_datadiff*0.05
    i_ymax = i_datamax + i_datadiff*0.05
    
    i_axis.vlines(events2019.loc[36:38].FormattedTime.values,i_ymin,i_ymax,linewidth=2,linestyle='--',color='#cc8904',label='Detection')
    i_axis.vlines(earthquakes.FormattedTime.values,i_ymin,i_ymax,linewidth=2,linestyle='-.',color='#5c04c7',label='Earthquake')

    ## Set ylimits:
    i_axis.set_ylim([i_ymin,i_ymax])
    
    ## Axis formatter:
    i_axis.yaxis.set_tick_params(labelsize=16)
    
    ## Grid:
    i_axis.grid()
    
    ## Legend:
    i_axis.legend(fontsize=16,loc='center right')

## x axis formatter etc at end        
## xlimits
i_axis.set_xlim([xmin,xmax])

i_axis.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d, %Hh'))
i_axis.xaxis.set_tick_params(labelsize=16,rotation=0)  

## label:
i_axis.set_xlabel('Month/Day/Hour in 2019',fontsize=17)
    

## Save:
multifig.savefig(main_multidataplot_path_pdf)
multifig.savefig(main_multidataplot_path_png)




# %%
## Make multi sensor plot with ADCP
multifig,multiaxes = plt.subplots(nrows=5,ncols=1,figsize=figure_dims,sharex=True)

dataframe_ordered = [oxygen_df,chlorophyll_df,temperature_df,turbidity_df]
column_ordered = ['oxygen','chlorophyll','seawatertemperature','turbidityntu']
ylabel_ordered = ['Oxygen (ml/L)','Chlorophyll (ug/L)','Sea Water \nTemperature (C)','Turbidity (NTU)']

## For each data type, plot:
for i_axis,i_dataframe,i_column,i_ylabel in zip(multiaxes[1:],dataframe_ordered,column_ordered,ylabel_ordered):
    ## get axis
    # i_axis = multiaxes[i_axis]
    ## scatter plot:
    i_axis.plot(i_dataframe.datetime,i_dataframe[i_column],linewidth=1.6,color='#056e6c')
    
    ## set label
    i_axis.set_ylabel(i_ylabel,fontsize=16)
    
    ## Add lines for events identified:
    i_datamin = i_dataframe.loc[((i_dataframe.datetime>xmin) & (i_dataframe.datetime<xmax))][i_column].min()
    i_datamax = i_dataframe.loc[((i_dataframe.datetime>xmin) & (i_dataframe.datetime<xmax))][i_column].max()
    i_datadiff = i_datamax - i_datamin
    
    i_ymin = i_datamin - i_datadiff*0.05
    i_ymax = i_datamax + i_datadiff*0.05
    
    i_axis.vlines(events2019.loc[36:38].FormattedTime.values,i_ymin,i_ymax,linewidth=2,linestyle='--',color='#cc8904',label='Detection')
    i_axis.vlines(earthquakes.FormattedTime.values,i_ymin,i_ymax,linewidth=2,linestyle='-.',color='#5c04c7',label='Earthquake')

    ## Set ylimits:
    i_axis.set_ylim([i_ymin,i_ymax])
    
    ## Axis formatter:
    i_axis.yaxis.set_tick_params(labelsize=16)
    
    ## Grid:
    i_axis.grid()
    
    ## Legend:
    i_axis.legend(fontsize=16,loc='center right')

## x axis formatter etc at end        
## xlimits
i_axis.set_xlim([xmin,xmax])

i_axis.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d, %Hh'))
i_axis.xaxis.set_tick_params(labelsize=16,rotation=0)  

## label:
i_axis.set_xlabel('Month/Day/Hour in 2019',fontsize=17)
    

## Save:
multifig.savefig(main_multidataplot_path_pdf)
multifig.savefig(main_multidataplot_path_png)

