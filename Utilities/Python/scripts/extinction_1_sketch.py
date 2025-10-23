"""
extinction_1_sketch.py
Converted from MATLAB script extinction_1_sketch.m to Python.

Makes figure for 'EXTINCTION 1' model for Tech Guide
"""

import os
import numpy as np
import fdsplotlib

plot_style = fdsplotlib.get_plot_style("fds")

# -----------------------------
# Paths (mirror extinction.py layout)
# -----------------------------
# Repository roots relative to this file
firemodels_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# Output directory consistent with MATLAB script intent
plot_dir = os.path.join(fds_dir, 'Manuals', 'FDS_Technical_Reference_Guide', 'FIGURES')
os.makedirs(plot_dir, exist_ok=True)

# -----------------------------
# Data (from MATLAB script)
# -----------------------------
# Temperature in °C
T1 = np.arange(0, 600 + 1, 1)       # 0:600
T2 = np.arange(600, 1500 + 1, 1)    # 600:1500

# Oxygen volume fraction pieces
X_O2_1 = 0.135 * (1.0 - (T1 - 20.0) / 1427.0)
X_O2_2 = 0.135 * (1.0 - (T2 - 20.0) / 1427.0)

# -----------------------------
# Make the figure (fdsplotlib style)
# -----------------------------
# Base figure with first segment (solid 'k-')
fig = fdsplotlib.plot_to_fig(
    x_data=T1,
    y_data=X_O2_1,
    marker_style='k-',
    linewidth=1,
    data_label=None,
    x_min=0.0, x_max=1500.0,
    y_min=0.0, y_max=0.2,
    xticks=[0, 500, 1000, 1500],
    yticks=np.linspace(0.0, 0.2, 5),
    x_label='Temperature (°C)',
    y_label='Oxygen Volume Fraction',
)

# Second segment (dashed 'k--')
fdsplotlib.plot_to_fig(
    x_data=T2,
    y_data=X_O2_2,
    marker_style='k--',
    linewidth=1,
    figure_handle=fig
)

# Vertical line at T = 600 °C from y=0 to y=0.08013
fdsplotlib.plot_to_fig(
    x_data=[600, 600],
    y_data=[0.0, 0.08013],
    marker_style='k-',
    linewidth=1,
    figure_handle=fig
)

# Axis label padding similar to extinction.py
ax = fig.axes[0]
ax.xaxis.set_tick_params(pad=10)
ax.yaxis.set_tick_params(pad=10)

# Text annotations (positions from MATLAB)
ax.text(150, 0.06, 'No Burn', fontsize=plot_style['Label_Font_Size'])
ax.text(800, 0.03, 'Burn', fontsize=plot_style['Label_Font_Size'])
ax.text(700, 0.15, 'Burn', fontsize=plot_style['Label_Font_Size'])

# Output path
out_pdf = os.path.join(plot_dir, 'extinction_1_sketch.pdf')

fig.savefig(out_pdf)
print(f"Saved figure to {out_pdf}")
