#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 15:13:18 2023

@author: vjs
"""
import numpy as np
import matplotlib.pyplot as plt
import onctools as ont
import pandas as pd    
from datetime import datetime, timedelta
import obspy as obs
import matplotlib.dates as mdates

# %%

##### PARAMETERS and paths

## Csv with list of templates, and flag if it is an event or not
template_list_path_base = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu'

template_list_paths_year = np.array([2018,2019,2020,2021,2022])

# Make list with template names, and master dir for miniseeds:
template_list_paths = []
template_masterdir_bases = []
for year in template_list_paths_year:
    ## get the list of paths for csvs:
    i_template_list_path = '%s%i/template_interpretations_BACAX%i.csv' % (template_list_path_base,year,year)
    template_list_paths.append(i_template_list_path)
    
    ## and list for master dirs:
    i_template_master_dir = '%s%i/' % (template_list_path_base,year)
    template_masterdir_bases.append(i_template_master_dir)
    

## Directory for figures
template_fig_dir = '/Users/vjs/turbidites/observational/figs/'




##### OTHER PARAMS ####

## PLOT params
plot_ncols = 5
## Plotsize:
plot_width = 16
subplot_height = 2.5 ## Width for each supbplot to get the total plot width

# %%
############################################################################## 

## For each year, read in the template list path, and append:
## Make empty dataframe:
template_list = pd.DataFrame()
for i_template_list_path, i_template_master_dir in zip(template_list_paths,template_masterdir_bases):
    i_template_list = pd.read_csv(i_template_list_path,skiprows=6,header=0)

    ## Convert the time columns to datetime:
    i_template_list['starttime'] = pd.to_datetime(i_template_list['time_min'], format='%Y%m%dT%H%M%S')
    i_template_list['endtime'] = pd.to_datetime(i_template_list['time_max'], format='%Y%m%dT%H%M%S')
    i_template_list['eventtime'] = pd.to_datetime(i_template_list['time_event'], format='%Y%m%dT%H%M%S')

    # Get the local paths for al template miniseeds
    print('iterating through templates...')

    # Make a list with the template full paths:
    i_template_full_paths = np.zeros_like(i_template_list['fname'].values)

    # Iterate through, get the turbidity template name, as well as the start time
    for index,row in i_template_list.iterrows():
        print('On template {}'.format(row.fname))
        
        ## Get the file base name:
        i_basename = i_template_list.loc[index]['fname'].split('/')[-1].split('.csv')[0]
        
        ## Get the path for the template:
        i_template_path = i_template_master_dir + 'csv2mseed/filtered/' + i_basename + '.mseed'
        
        i_template_full_paths[index] = i_template_path
        
        print(i_template_path)
        
    ## Add to index's dataframe:
    i_template_list['local_fname'] = i_template_full_paths

    ## Append this entire dataframe to the original:
    template_list = pd.concat([template_list,i_template_list],ignore_index=True)



# %%
## Make a plot with all the templates on it

## Get the number of rows that will be needed:
templates_plotDF = template_list.loc[template_list['TurbidityEvent_woADCP'] == True].reset_index(drop='True')
num_templates = len(templates_plotDF)

## How many rows with ncols?
plot_nrows = np.ceil(num_templates/plot_ncols).astype('int')


## Initiate plot:
template_fig, template_axes = plt.subplots(nrows=plot_nrows, ncols=plot_ncols)

if len(templates_plotDF) > plot_ncols:
    template_fig.set_figheight((subplot_height+0.1)*plot_nrows)
else:
    template_fig.set_figheight(subplot_height+0.1)

template_fig.set_figwidth(plot_width)

## Loop through templates and plot:
## Set a template counter:
template_counter = 0
for i_row in range(plot_nrows):
    for i_col in range(plot_ncols):
        ## If you're not past the total number of templates:
        if template_counter < len(templates_plotDF):
            ## Get the axes:
            if i_row > 0:
                i_axes = template_axes[i_row,i_col]
            elif i_row == 0 and len(templates_plotDF) <= plot_ncols:
                i_axes = template_axes[i_col]
            else:
                i_axes = template_axes[i_row,i_col]
            
            ## Get the template data:
            i_turbidity = obs.read(templates_plotDF.loc[template_counter].local_fname)
                        
            
            ## Plot:
            i_axes.plot(i_turbidity[0].times('matplotlib'),i_turbidity[0].data,color='#145246',linewidth=1.5)
                
            ## Make sure x axis is UTC Datetime:
            i_axes.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            
            # Rotates and right-aligns the x labels so they don't crowd each other.
            for label in i_axes.get_xticklabels(which='major'):
                label.set(rotation=30, horizontalalignment='right', fontsize=11)
            
            ## Add a grid:
            i_axes.grid()
                            
            ## Add a title, x axis label, y axis label
            i_title = templates_plotDF.time_event[template_counter]
            i_axes.set_title(i_title)
            
            ## If it's the last row, add an x label:
            if i_row == plot_nrows-1:
                i_axes.set_xlabel('Time HH:MM',fontsize=12)
            
            ## IF it's the second to last row but no panels underneath:
            modulo_of_last_column = len(templates_plotDF) % plot_ncols
            
            if ((i_row == plot_nrows-2) & (i_col >= modulo_of_last_column)):
                i_axes.set_xlabel('Time HH:MM',fontsize=12)
            
            ## IF it's on the left hand side, add a y label:
            if i_col == 0:
                i_axes.set_ylabel('Turbidity (NTU)',fontsize=12)
            
                            
            ## Add to counter:
            template_counter += 1
            
            
        ## Otherwise remove this axis:
        else:
            if i_row > 0:
                template_fig.delaxes(template_axes[i_row][i_col])
            else:
                template_fig.delaxes(template_axes[i_col])
    
## Adjust:
if len(templates_plotDF) > plot_ncols:
    template_fig.subplots_adjust(wspace=0.3,hspace=0.38,left=0.07,right=0.93,top=0.93,bottom=0.07)
else:
    template_fig.subplots_adjust(left=0.07,right=0.93,top=0.89,bottom=0.28)


## Save figure:
template_fig.savefig(template_fig_dir+'pdf/all_templates_woADCP.pdf',transparent=True)
template_fig.savefig(template_fig_dir+'png/all_templates_woADCP.png',transparent=True)
    
