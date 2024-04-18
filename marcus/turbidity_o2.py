#!/usr/bin/env python
# coding: utf-8

# ## Multivariable Plots
# Makes Multivariable plots as a function of time.
# 
# Input: 
# - Turbidity file, with temperature
# - Oxygen file
# 
# Output:
# - Creates, but doesn't save multivariable plot
# 
# Three subplots (turbidity, O2, and temperature).

# In[76]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.dates as mdates


# In[103]:


data= pd.read_csv('/Users/jabba7hebutt/Downloads/search25234278/BarkleyCanyon_BarkleyCanyonHead_TurbidityMeter_20181022T000000Z_20181023T235500Z-NaN_clean.csv',skiprows=52)
data1 = pd.read_csv('/Users/jabba7hebutt/Downloads/search25234281/BarkleyCanyon_BarkleyCanyonHead_OxygenSensor_20181022T000000Z_20181024T000000Z-NaN_clean_avg1hour.csv',skiprows=52)


# In[104]:


time = np.array(data['#"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"'])
time1 = np.array(data1['#"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"'])
turb = np.array(data[' "Turbidity (FTU)"'])
#chloro = np.array(data[' "Chlorophyll (ug/l)"'])
o2 = np.array(data1[' "Oxygen Concentration Corrected (ml/l)"'])
temp = np.array(data[' "Temperature (C)"'])


# In[113]:


sns.set()
fig, axes = plt.subplots(3,sharex=True)
fig.set_size_inches(20, 9)

p0 = axes[0].plot(time,turb)
axes[0].set_ylabel('Turbidity (FTU)',fontsize=18)

p1 = axes[1].plot(time,temp)
axes[1].set_ylabel('Temp (C)',fontsize=18)

p2 = axes[2].plot(time1,o2)
axes[2].set_ylabel('O2 conc. (ml/l)',fontsize=18)

axes[2].xaxis.set_major_locator(mdates.DayLocator(interval=40))

plt.gcf().autofmt_xdate()
plt.show()


# In[ ]:




