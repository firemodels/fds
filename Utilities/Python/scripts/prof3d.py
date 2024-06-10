"""
McDermott
6-10-2024

Animate a 3D scatterplot of FDS in-depth profile (CHID_prof_n.csv)

Usage: Copy this script to your working directory, change the filename.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Close all previously opened figures
plt.close('all')

tmpa = 20.

filename = 'pine_21O2_40_1C_cat_prof_1.csv'

# read header information

with open(filename) as f:
    first_line = f.readline().strip('\n')

XYZ=first_line.split(",")[1:4]
X = float(XYZ[0]) #; print(X)
Y = float(XYZ[1]) #; print(Y)
Z = float(XYZ[2]) #; print(Z)

df = pd.read_csv(filename,skiprows=2,header=None)
df.fillna(tmpa, inplace=True)

t = df.values[:,0]
time_diff = np.diff(t)  # Calculate time difference between consecutive frames

n_points = len(t)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the scatter plot (initially empty)
scatter = ax.scatter([], [], [], c=[], cmap=plt.hot(), vmin=20, vmax=500)
cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=5)
cbar.set_label('Scalar Values')

def update(frame):
    n = int(df.values[frame,1])
    x = df.values[frame,2:2+n]
    xs = X + x # in-depth positions could update during simulation
    ys = Y + np.zeros(len(x))
    zs = Z + np.zeros(len(x))
    scalar_values = df.values[frame,2+n:2+2*n]

    # Clear the previous plot
    ax.cla()
    # Plot the scatter plot for the current frame
    scatter = ax.scatter(xs, ys, zs, c=scalar_values, cmap=plt.hot(), marker='o')

    # Set appropriate labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # Display current time
    ax.text(0,0,0.02,'Time: {:.2f}'.format(t[frame]),size=16, zorder=1,color='black')

# Create the animation
interval = 200 # milliseconds
ani = animation.FuncAnimation(fig, update, frames=n_points, interval=interval, blit=False, repeat=False)


plt.show()



