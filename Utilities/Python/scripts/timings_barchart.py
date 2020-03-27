#!/usr/bin/env python3
# McDermott
# March 2020

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')

B = pd.read_csv('methane_dx_p3125cm_cpu.csv', sep=' *, *', engine='python', comment='#')
G = pd.read_csv('propane_dx_p3125cm_cpu.csv', sep=' *, *', engine='python', comment='#')

labels = ['MAIN', 'DIVG', 'MASS', 'VELO', 'PRES', 'WALL', 'DUMP', 'RADI', 'FIRE', 'COMM', 'Total T_USED (s)']

h_burn = B[labels].values[10]/3600. # look at the 10th output
h_gord = G[labels].values[10]/3600. # look at the 10th output

x = np.arange(len(labels))
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, h_burn, width, label='burn')
rects2 = ax.bar(x + width/2, h_gord, width, label='gordon')

ax.set_ylabel('time (hrs)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
plt.xticks(rotation=45)

plt.tight_layout()
plt.show()
plt.savefig('timings.pdf')