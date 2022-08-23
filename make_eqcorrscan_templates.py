#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 13:15:52 2022

@author: vjs
"""
import logging

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.core import template_gen

import obspy as obs

import numpy as np
import pandas as pd

import onc_download_module as odm

from glob import glob

import matplotlib.pyplot as plt

#################################
# %% PARAMTERS/PATHS

## Template csv directory:
csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/'
    
## Template output mseed directory:
mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/mseed/'

## Template figure directory:
mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/plots/mseed/'

## DAta type/column name for htese templates:
data_type = 'turbntu'

## Unit of data:
data_unit = 'NTU'

## dta network:
data_network = 'ONC'
station_name = 'BACAX'

## Figure size:
plotsize = (8,4)

    
# %%
## Loop through the files in teh template directory:
all_csv_paths = glob(csv_directory + '*.csv')

## Loop through them:
for i_template in range(len(all_csv_paths)):
    i_template_path = all_csv_paths[i_template]
    
    print('working on %s \n \n ' % i_template_path)
    
    i_template_df = pd.read_csv(i_template_path,names=['datetime_string',data_type])

    ## Convert datetime string to datetime:
    i_template_df['datetime'] = i_template_df.datetime_string.astype('datetime64[ms]')

    ## convert to a stream object
    i_stream = odm.convert_df2stream(i_template_df,data_type,datetime_name='datetime',network=data_network,station=station_name,location='',channel=data_type)

    ## Save it to file:
    i_mseed_base = i_template_path.split('/')[-1].split('.csv')[0]
    i_mseed_path = mseed_directory + i_mseed_base + '.mseed'
    i_stream.write(i_mseed_path,format='MSEED')

    ## Plot:
    i_plot_path = mseed_figure_directory + i_mseed_base + '.png'
    
    ## INitiate plot and other commands before plotting:
    i_fig, i_ax = plt.subplots(figsize=plotsize)
    i_ax.set_ylabel(data_unit)
    i_ax.xaxis.set_ticklabels([])
    i_ax.yaxis.set_ticklabels([])
    
    i_stream.plot(fig=i_fig)

    
    print('saving to %s \n \n \n' % i_plot_path)
    ## Save plot:
    i_fig.savefig(i_plot_path)

    ## Close plot:
    plt.close('all')



# %% OLD TEST
# ## Read in a template from Debi:
# test_template_path = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/reformated_BACAX_ntu_2019_template_notfiltered_2_15_1.1254_1_0_0.5_20191224T035330.csv'

# test_template = pd.read_csv(test_template_path,names=['datetime_string','turbidityntu'])

# ## Convert datetime string to datetime:
# test_template['datetime'] = test_template.datetime_string.astype('datetime64[ms]')

# ## Make into an obspy stream - first into a trace:
# # Fill header attributes
# stats = {'network': 'ONC', 'station': 'BACAX', 'location': '',
#          'channel': 'turbidityntu', 'npts': len(test_template.turbidityntu), 'sampling_rate': 0.1,
#          'mseed': {'dataquality': 'D'}}



# dataframe = test_template
# data_name = 'turbidityntu'
# datetime_name = 'datetime'
# network = 'ONC'
# station = 'BACAX'
# channel = 'turbidityntu'


# test_indices = np.r_[np.arange(25),np.arange(30,70),np.arange(75,86),np.arange(100,122),np.arange(130,161)]
# new_test_template = test_template.loc[test_indices] #.reset_index(dropna=True)

# test_stream = odm.convert_df2stream(new_test_template,data_name,datetime_name='datetime',network='ONC',station='BACAX',location='',channel='turbidityntu')




