#!/usr/bin/env python3
# McDermott
# Oct 2024

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')

wrk_dir = '../../../Verification/Fires/'

M = pd.read_csv(wrk_dir+'simple_test_hrr.csv', sep=' *, *', engine='python', comment='#', header=1)

time   = M['Time'].values
hrr    = M['HRR'].values
q_radi = M['Q_RADI'].values

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('hrr [kW]', color=color)
ax1.plot(time, hrr, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('q_radi [kW]', color=color)  # we already handled the x-label with ax1
ax2.plot(time, q_radi, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# ax.legend()

# plt.savefig('plot_2yaxes_example.pdf')