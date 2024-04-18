#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 10:21:22 2022

@author: vjs
"""

from onc.onc import ONC
import requests
import re
import tkinter as tk
from appdirs import user_data_dir
import os
#root = tk.Tk()  
import pandas as pd
import requests
import json
import os
from contextlib import closing
import errno
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

# %% Paths
barkley_codes_path = '/Users/vjs/turbidites/observational/data/onc_codes/locationcodes_barkley.csv'
clayoquot_codes_path = '/Users/vjs/turbidites/observational/data/onc_codes/locationcodes_clayoquot.csv'


turbidity_codes_path = '/Users/vjs/turbidites/observational/data/onc_codes/propertycodes_turbidity.csv'



## Test download path:
temp_dL_dir = '/Users/vjs/turbidites/observational/data/temp'

## codes etc. to search for
propertyCodes = ['turbidityftu','turbidityntu','seawatertemperature','oxygen','pressure','chlorophyll']
locationCodes = ['BACAX']

# %% Token function
def Token():
    global tokenpath
    tokenpath = user_data_dir("ONC-Token", "ONC")
    if os.path.exists(tokenpath + r"\token.txt"):
        print("Token file exists.")
        f = open(tokenpath + r"\token.txt", "r")
        token = f.read()
        f.close()
        return token
    

# %%  Get location codes for Barkley canyon %% #
token = Token()
onc = ONC(token)
 

# %% Download location codes 

## Barkley Canyon:
barkley_filters = {'locationName': ['Bark']}
barkley_locationcodes_result = onc.getLocations(barkley_filters)
 
onc.print(barkley_locationcodes_result)

## Make a dataframe:
barkleydf = pd.DataFrame(barkley_locationcodes_result)
## SAve to file
barkleydf.to_csv(barkley_codes_path)


## Clayoquot Slope:
clayoquot_filters = {'locationName': ['Clayoquot','Bullseye','Bubbly','Gastown','ODP']}
clayoquot_locationcodes_result = onc.getLocations(clayoquot_filters)
 
onc.print(clayoquot_locationcodes_result)

## Make a dataframe:
clayoquotdf = pd.DataFrame(clayoquot_locationcodes_result)
## SAve to file
clayoquotdf.to_csv(clayoquot_codes_path)



# %% Get property codes for turbidity, ctd, oxygen

## Make a filter that will search for property codes fo rthings that have 
##   'turbidity' in the description
turbidity_filters = {'description':'turbidity'}

## search for property codes
turbidity_propCodes_result = onc.getProperties(turbidity_filters)

## Make a dataframe:
turbiditydf = pd.DataFrame(turbidity_propCodes_result)
## Save:
turbiditydf.to_csv(turbidity_codes_path)


# %% Try finding data at a given location code, and proprty code

url = 'https://data.oceannetworks.ca/api/dataProducts'
parameters = {'method':'get',
            'token':token, # replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in.
            'locationCode':'BACAX',
            'propertyCode':'turbidityftu'}
  
response = requests.get(url,params=parameters)
  
if (response.ok):
    dataProducts = json.loads(str(response.content,'utf-8')) # convert the json response to an object
    for dataProduct in dataProducts:
        print(dataProduct)
else:
    if(response.status_code == 400):
        error = json.loads(str(response.content,'utf-8'))
        print(error) # json response contains a list of errors, with an errorMessage and parameter
    else:
        print ('Error {} - {}'.format(response.status_code,response.reason))


# %% Ask for a data product delivery, given search paramters. Makes a response which
##   is then Ran in a Data Run Request...
url = 'https://data.oceannetworks.ca/api/dataProductDelivery'
parameters = {'method':'request',
            'token':token,          # replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in.
           'locationCode':'BACAX',             # Barkley Canyon / Axis (POD 1)
            'propertyCode':propertyCodes[1],    # 150 kHz Acoustic Doppler Current Profiler
            'dataProductCode':'TSSD',           # Time Series Scalar Data
            'extension':'csv',                  # Comma Separated spreadsheet file
            'dateFrom':'2019-12-24T00:00:00.000Z',  # The datetime of the first data point (From Date)
            'dateTo':'2019-12-26T00:00:00.000Z',    # The datetime of the last data point (To Date)
            'dpo_qualityControl':1,             # The Quality Control data product option - See https://wiki.oceannetworks.ca/display/DP/1
            'dpo_resample':'average',              # The Resampling data product option - See https://wiki.oceannetworks.ca/display/DP/1
            'dpo_average':60,                   # Resampling average, 60=1 minutes
            'dpo_dataGaps':0}                   # The Data Gaps data product option - See https://wiki.oceannetworks.ca/display/DP/1
 
response = requests.get(url,params=parameters)
  
if (response.ok):
    requestInfo = json.loads(str(response.content,'utf-8')) # convert the json response to an object
     
    print('Request Id: {}'.format(requestInfo['dpRequestId']))      # Print the Request Id
     
    if ('numFiles' in requestInfo.keys()):
        print('File Count: {}'.format(requestInfo['numFiles']))     # Print the Estimated File Size
  
    if ('fileSize' in requestInfo.keys()):
        print('File Size: {}'.format(requestInfo['fileSize']))      # Print the Estimated File Size
     
    if 'downloadTimes' in requestInfo.keys():
        print('Estimated download time:')
        for e in sorted(requestInfo['downloadTimes'].items(),key=lambda t: t[1]):
            print('  {} - {} sec'.format(e[0],'{:0.2f}'.format(e[1])))
 
 
    if 'estimatedFileSize' in requestInfo.keys():
        print('Estimated File Size: {}'.format(requestInfo['estimatedFileSize']))
                 
    if 'estimatedProcessingTime' in requestInfo.keys():
        print('Estimated Processing Time: {}'.format(requestInfo['estimatedProcessingTime']))
  
else:
    if(response.status_code == 400):
        error = json.loads(str(response.content,'utf-8'))
        print(error) # json response contains a list of errors, with an errorMessage and parameter
    else:
        print ('Error {} - {}'.format(response.status_code,response.reason))
        
# %%
## Run the data product request...
url = 'https://data.oceannetworks.ca/api/dataProductDelivery'      
parameters = {'method':'run',
            'token':token,              # replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in.
            'dpRequestId':requestInfo['dpRequestId']}     # replace YOUR_REQUEST_ID_HERE with a requestId number returned from the request method
response = requests.get(url,params=parameters)
  
if (response.ok):
    r = json.loads(str(response.content,'utf-8')) # convert the json response to an object
  
    for runId in [run['dpRunId'] for run in r]:
        print('Run Id: {}'.format(runId))       # Print each of the Run Ids
 
else:
    if(response.status_code == 400):
        error = json.loads(str(response.content,'utf-8'))
        print(error) # json response contains a list of errors, with an errorMessage and parameter
    else:
        print ('Error {} - {}'.format(response.status_code,response.reason))
        
# %% Download the data. Saves to outPath
url = 'https://data.oceannetworks.ca/api/dataProductDelivery'
parameters = {'method':'download',
            'token':token,   # replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in..
            'dpRunId':r[0]['dpRunId'],       # replace YOUR_RUN_ID with the dpRunId returned from the 'run' method.

            'index':1}                   # for run requests that contain more than one file, change the index number to the index of the file you would like to download.
                                           # If the index number does not exist an HTTP 410 and a message will be returned.
 
 
outPath=temp_dL_dir                        # replace with the file location you would like the file to be downloaded to.

status_code_check = False 

while status_code_check == False:
    print('trying download')
    with closing(requests.get(url,params=parameters,stream=True)) as streamResponse:
        if streamResponse.status_code == 200: #OK
            ## Change code to false so it doesn't try to download again
            status_code_check = True
            print('it worked')
            if 'Content-Disposition' in streamResponse.headers.keys():
                content = streamResponse.headers['Content-Disposition']
                filename = content.split('filename=')[1]
            else:
                print('Error: Invalid Header')
                streamResponse.close()
                sys.exit(-1)
             
            filePath = '{}/{}'.format(outPath,filename)
            try:
                if (not os.path.isfile(filePath)):
                    #Create the directory structure if it doesn't already exist
                    try:
                        os.makedirs(outPath)
                    except OSError as exc:
                        if exc.errno == errno.EEXIST and os.path.isdir(outPath):
                            pass
                        else:
                            raise
                    print ("Downloading '{}'".format(filename))
     
                    with open(filePath,'wb') as handle:
                        try:
                            for block in streamResponse.iter_content(1024):
                                handle.write(block)
                        except KeyboardInterrupt:
                            print('Process interupted: Deleting {}'.format(filePath))
                            handle.close()
                            streamResponse.close()
                            os.remove(filePath)
                            sys.exit(-1)
                else:
                    print ("  Skipping '{}': File Already Exists".format(filename))
            except:
                msg = 'Error streaming response.'
                print(msg)
        else:
            if(streamResponse.status_code in [202,204,400,404,410]):
                payload = json.loads(str(streamResponse.content,'utf-8'))
                if len(payload) >= 1:
                    msg = payload['message']
                    print('HTTP {} - {}: {}'.format(streamResponse.status_code,streamResponse.reason,msg))
            else:
                print ('Error {} - {}'.format(streamResponse.status_code,streamResponse.reason))
     
    streamResponse.close()

# %%
## Try to plot...

temppath = '/Users/vjs/turbidites/observational/data/temp.csv/BarkleyCanyon_BarkleyCanyonAxis_variables_TurbidityNTU_20191224T000000Z_20191226T000000Z-clean_avg1minute.csv'
testturbdata = pd.read_csv(temppath,skiprows=54,names=['datetime','turbidity','qcflag','turbiditycount'])

plt.plot(testturbdata.datetime,testturbdata.turbidity,linewidth=1,color='blue')
plt.show()

plt.xlabel('Date Time')
plt.ylabel('Turbidity...')
plt.savefig(temp_dL_dir + '/testdL_BACAX_Dec2019.png')
