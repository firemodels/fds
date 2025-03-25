"""
McDermott
6-10-2024

Animate a 3D scatterplot of FDS in-depth profile (CHID_prof_n.csv)

Usage: Copy this script to your working directory and add your
       profiles to the filenames list.

       Use option --with_slider to control the animation with a time slider
       Use option --save_animation to save the animation (no slider) to a movie file

Example:
$ python prof3d.py --with_slider
"""

import sys
import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--with_slider', action='store_true', help='Control animation with a time slider')
parser.add_argument('--save_animation', action='store_true', help='Save animation')

args = parser.parse_args()

# if args.with_slider:
#     print("Flag is present")
# else:
#     print("Flag is not present")

# sys.exit()

# Close all previously opened figures
plt.close('all')

scalar_min = 0.
scalar_max = 120.

filenames = ['../Current_Results/pine_21O2_40_1C_cat_prof_4.csv']

# create lists to store information about each profile
IOR = []
X = []
Y = []
Z = []
df = {}

for i, filename in enumerate(filenames):

    data = []
    max_cols = 0

    # read header information

    with open(filename,'r') as f:
        # Skip the first 1 lines
        for j in range(1):
            next(f)
        first_line = f.readline().strip('\n')

        header=first_line.split(",")[1:5]
        IOR.append(int(header[0])) ; print(IOR)
        X.append(float(header[1])) ; print(X)
        Y.append(float(header[2])) ; print(Y)
        Z.append(float(header[3])) ; print(Z)
        next(f)

        # Read lines one at a time
        while True:
            line = f.readline()
            if not line:  # End of file
                break
            row = line.strip().split(',')  # Adjust delimiter if needed

            # Convert each element to float, handling errors gracefully
            try:
                row = [float(value) if value else None for value in row]
            except ValueError:
                # If a value cannot be converted, keep it as None
                row = [float(value) if value.replace('.', '', 1).isdigit() else None for value in row]

            data.append(row)
            max_cols = max(max_cols, len(row))  # Track the maximum number of columns

    # Normalize rows to have the same number of columns
    data = [row + [None] * (max_cols - len(row)) for row in data]

    # Convert to a Pandas DataFrame
    df[i] = pd.DataFrame(data)

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

if args.with_slider:
    # Create a slider for controlling time
    axslider = plt.axes([0.1, 0.02, 0.8, 0.03])
    time_slider = Slider(axslider, 'Time', t[0], t[-1], valinit=t[0])

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

if args.with_slider:
    # Function to update the plot when the slider is moved
    def update_slider(val):
        time = time_slider.val
        index = np.abs(t - time).argmin()  # Find the index closest to the current time
        update(index)
        fig.canvas.draw_idle()  # Redraw the plot

    # Connect the slider to the update_slider function
    time_slider.on_changed(update_slider)
else:
    # Create the animation
    interval = 200 # milliseconds
    ani = animation.FuncAnimation(fig, update, frames=n_points, interval=interval, blit=False, repeat=False)

    if args.save_animation:
        # Save the animation as an mp4 file
        ani.save('animation.mp4', writer='ffmpeg')

plt.show()



