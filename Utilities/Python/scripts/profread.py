"""
McDermott
5-24-2024

Read and plot FDS in-depth profile (CHID_prof_n.csv)

Usage: Copy this script to your working directory, change the filename in pd.read_csv, and choose a time index
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('pine_21O2_40_1C_cat_prof_1.csv',skiprows=3)

t = df.values[:,0]
n = int(df.values[1,1])

t_index = -1 # last value

# # or specify a time (comment this block to retain last value)
# specific_time = 600.
# difference_array = np.absolute(t-specific_time) # calculate the difference array
# t_index = difference_array.argmin() # find the index of minimum element from the array (why this cannot be done in one step is beyond me)

x = df.values[t_index,2:2+n]
y = df.values[t_index,2+n:3+2*n]

fig, ax = plt.subplots()

ax.plot(x,y,marker='o',linestyle='-',linewidth=1)

# set your axis labels
plt.xlabel('Position (m)')
plt.ylabel('Temperature ($^\circ$C)')

plt.show()

# plt.savefig('profile.pdf')