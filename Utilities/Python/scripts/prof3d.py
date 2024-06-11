"""
McDermott
6-10-2024

Animate a 3D scatterplot of FDS in-depth profile (CHID_prof_n.csv)

Usage: Copy this script to your working directory, change the filename.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Close all previously opened figures
plt.close('all')

tmpa = 20.
scalar_min = 20.
scalar_max = 60.

filenames = ['pine_21O2_40_1C_cat_prof_1.csv',
             'pine_21O2_40_1C_cat_prof_5.csv',
             'pine_21O2_40_1C_cat_prof_9.csv',
             'pine_21O2_40_1C_cat_prof_13.csv']

# create lists to store information about each profile
IOR = []
X = []
Y = []
Z = []
df = []

for i in range(len(filenames)):

    # read header information

    with open(filenames[i]) as f:
        # Skip the first 1 lines
        for j in range(1):
            next(f)
        first_line = f.readline().strip('\n')

    header=first_line.split(",")[1:5]
    IOR.append(int(header[0])) #; print(IOR)
    X.append(float(header[1])) #; print(X)
    Y.append(float(header[2])) #; print(Y)
    Z.append(float(header[3])) #; print(Z)

    # sys.exit()

    df.append(pd.read_csv(filenames[i],skiprows=3,header=None))
    df[i].fillna(tmpa, inplace=True)

# sys.exit()

# take time values from first profile

t = df[0].values[:,0]
time_diff = np.diff(t)  # Calculate time difference between consecutive frames

n_points = len(t)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
colormap = plt.jet() # see https://matplotlib.org/stable/gallery/color/colormap_reference.html

# Plot the scatter plot (initially empty)
scatter = ax.scatter([], [], [], c=[], cmap=colormap, vmin=scalar_min, vmax=scalar_max)
cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5)
cbar.set_label('Scalar Values')

def update(frame):
    # Clear the previous plot
    ax.cla()

    for i in range(len(filenames)):
        n = int(df[i].values[frame,1])
        x = df[i].values[frame,2:2+n] # in-depth positions
        xs = X[i] + np.zeros(len(x))
        ys = Y[i] + np.zeros(len(x))
        zs = Z[i] + np.zeros(len(x))
        match IOR[i]:
            case -1: xs = X[i] + x
            case  1: xs = X[i] - x
            case -2: ys = Y[i] + x
            case  2: ys = Y[i] - x
            case -3: zs = Z[i] + x
            case  3: zs = Z[i] - x

        scalar_values = df[i].values[frame,2+n:2+2*n]

        # Plot the scatter plot for the current frame
        scatter = ax.scatter(xs, ys, zs, c=scalar_values, cmap=colormap, marker='.', s=100, vmin=scalar_min, vmax=scalar_max, alpha=1, depthshade=False)

    # Set appropriate labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # Display current time
    ax.text2D(0.05, 0.95, 'Time: {:.2f}'.format(t[frame]), size=16, zorder=1, transform=ax.transAxes)

# Create the animation
interval = 200 # milliseconds
ani = animation.FuncAnimation(fig, update, frames=n_points, interval=interval, blit=False, repeat=False)

plt.show()



