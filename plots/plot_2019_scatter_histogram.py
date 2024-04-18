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

## Output figure path:
jointplot_path_pdf = '/Users/vjs/turbidites/observational/figs/pdf/detections2019jointplot.pdf'
jointplot_path_png = '/Users/vjs/turbidites/observational/figs/png/detections2019jointplot.png'


### Read it in:
# %%

## Read in turbidity events identified by STA/LTA:
eventsall = pd.read_csv(datapath,names=['Time','Amp','Review'],skiprows=1)
## Get datetime column:
eventsall['FormattedTime'] = pd.to_datetime(eventsall.Time.astype('datetime64[ms]'))
## And take only vetted ones:
eventdata = eventsall.loc[eventsall.Review == 1]

## Dictionary for histogram plots for color etc.
histformat_dict = dict(bins=12,color='#16627d')

## make plot
eventplot = sns.jointplot(data=eventdata,x='FormattedTime',y='Amp',alpha=0.5,s=500,facecolor='#16627d',edgecolor='black',marginal_kws=histformat_dict,height=12)

## Labels:
eventplot.ax_joint.set_xlabel('Time of Detection, 2019',fontsize=28)
eventplot.ax_joint.set_ylabel('Turbidity Event Amplitude (NTU)',fontsize=28)
eventplot.fig.subplots_adjust(bottom=0.15,left=0.1)

## Ticks and params
eventplot.ax_joint.xaxis.set_major_formatter(mdates.DateFormatter('%B'))
eventplot.ax_joint.xaxis.set_tick_params(labelsize=19,rotation=0)
eventplot.ax_joint.yaxis.set_tick_params(labelsize=24)

## ADd grid:
eventplot.ax_joint.grid()

## Axis limits:
eventplot.ax_joint.set_xlim([datetime(2018,12,29),datetime(2019,12,31)])

## Save plot:
eventplot.savefig(jointplot_path_pdf,transparent=True)
eventplot.savefig(jointplot_path_png,transparent=True)


