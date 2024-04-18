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



#############
## parameters and paths

datapath = '/Users/vjs/turbidites/observational/data/templates/events_turbidity_BACAX2019.csv'

template_path_base = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu'

year_list = [2018,2019,2020,2021]


## Output figure path:
jointplot_path_pdf = '/Users/vjs/turbidites/observational/figs/pdf/detections_2018_2021jointplot.pdf'
jointplot_path_png = '/Users/vjs/turbidites/observational/figs/png/detections_2018_2021jointplot.png'


# %%
## GEt list of paths to read in:
template_path_list = []
for i_year in year_list:
    i_path = template_path_base + np.str(i_year) + '/template_interpretations_BACAX' + np.str(i_year) + '.csv'    
    template_path_list.append(i_path)
    
# %%
### Read it in:

for i_year_index,i_path in enumerate(template_path_list):
    if i_year_index == 0:
        template_list = pd.read_csv(i_path,skiprows=6,header=0)
    else:
        i_template_list = pd.read_csv(i_path,skiprows=6,header=0)
        
        ## Add to the template dataframe:
        template_list = pd.concat([template_list,i_template_list],ignore_index=True)

## Convert the time columns to datetime:
template_list['starttime'] = pd.to_datetime(template_list['time_min'], format='%Y%m%dT%H%M%S')
template_list['endtime'] = pd.to_datetime(template_list['time_max'], format='%Y%m%dT%H%M%S')
template_list['eventtime'] = pd.to_datetime(template_list['time_event'], format='%Y%m%dT%H%M%S')

    
## Take only "true" ones:
eventdata = template_list.loc[template_list.TurbidityEvent_woADCP == True].reset_index(drop='True')

## Try grouping by month for plot:
# create a representation of the month with strfmt
eventdata['month_of_event'] = eventdata['eventtime'].map(lambda dt: dt.strftime('%m')).astype('int')



# %% Plotting
## Dictionary for histogram plots for color etc.
histformat_dict = dict(bins=12,color='#16627d')

## make plot
# eventplot = sns.jointplot(data=eventdata,x='FormattedTime',y='Amp',alpha=0.5,s=500,facecolor='#16627d',edgecolor='black',marginal_kws=histformat_dict,height=12)
eventplot = sns.jointplot(data=eventdata,x='month_of_event',y='amp_event',alpha=0.5,s=500,facecolor='#16627d',edgecolor='black',marginal_kws=histformat_dict,height=12)


## Labels:
eventplot.ax_joint.set_xlabel('Month of Detection',fontsize=28)
eventplot.ax_joint.set_ylabel('Turbidity Event Amplitude (NTU)',fontsize=28)
eventplot.fig.subplots_adjust(bottom=0.15,left=0.1)

## Ticks and params
#eventplot.ax_joint.xaxis.set_major_formatter(mdates.DateFormatter('%B'))
eventplot.ax_joint.xaxis.set_tick_params(labelsize=19,rotation=0)
eventplot.ax_joint.yaxis.set_tick_params(labelsize=24)

## ADd grid:
eventplot.ax_joint.grid()

## Axis limits:
eventplot.ax_joint.set_xlim([0,13])


## Save plot:
eventplot.savefig(jointplot_path_pdf,transparent=True)
eventplot.savefig(jointplot_path_png,transparent=True)


