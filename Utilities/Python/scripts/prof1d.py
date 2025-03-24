"""
McDermott
12-20-2024

Animate a 1D plot of FDS in-depth profile (CHID_prof_n.csv)

Usage: Copy this script to your working directory and add your
       profiles to the filenames list.

       Use option --with_slider to control the animation with a time slider
       Use option --save_animation to save the animation (no slider) to a movie file

Example:

$ python prof1d.py --with_slider

"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
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

scalar_min = [20  ,  0,   0,   0,  0, 0]
scalar_max = [1000, 20, 400, 120, 10, 1]
scalar_lab = ['T','H2O','PINE','CHAR','ASH','O2']

filenames = ['../Test/pine_21O2_40_1C_cat_prof_1.csv',
             '../Test/pine_21O2_40_1C_cat_prof_2.csv',
             '../Test/pine_21O2_40_1C_cat_prof_3.csv',
             '../Test/pine_21O2_40_1C_cat_prof_4.csv',
             '../Test/pine_21O2_40_1C_cat_prof_5.csv',
             '../Test/pine_21O2_40_1C_cat_prof_6.csv']

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
        IOR.append(int(header[0])) #; print(IOR)
        X.append(float(header[1])) #; print(X)
        Y.append(float(header[2])) #; print(Y)
        Z.append(float(header[3])) #; print(Z)
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

fig, ax = plt.subplots(3, 2, figsize=(12, 8))
ax = ax.flatten()  # Flatten the array for easier indexing

# Initialize the figure and axis
# fig, ax = plt.subplots()
for i in range(len(filenames)):
    line, = ax[i].plot([], [], marker='o', label="Scalar vs X")

time_text = fig.text(0.5, 0.95, "", ha="center", va="center", fontsize=14, weight="bold")

if args.with_slider:
    # Create a slider for controlling time
    axslider = plt.axes([0.1, 0.02, 0.8, 0.03])
    time_slider = Slider(axslider, 'Time', t[0], t[-1], valinit=t[0])

def update(frame):

    for i in range(len(filenames)):
        n = int(df[i].values[frame,1])
        x = df[i].values[frame,2:2+n]*100. # in-depth positions (cm)
        scalar_values = df[i].values[frame,2+n:2+2*n]

        # Clear the previous plot
        ax[i].cla()

        # Plot the scatter plot for the current frame
        ax[i].plot(x, scalar_values, marker='o')
        if scalar_lab[i]=='O2' and np.max(scalar_values)>scalar_min[i]:
            ax[i].set_ylim(scalar_min[i], np.max(scalar_values)*1.1)
        else:
            ax[i].set_ylim(scalar_min[i], scalar_max[i])

        # Set appropriate labels
        ax[i].set_ylabel(scalar_lab[i])

        # Only show X-axis labels and ticks on the bottom row
        if i >= len(filenames) - 2:  # Adjust based on layout (3 rows, 2 columns)
            ax[i].set_xlabel('x (cm)')
        else:
            ax[i].set_xticklabels([])  # Remove labels

    # Update the shared time text
    time_text.set_text('Time: {:.2f}'.format(t[frame]))

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
    ani = FuncAnimation(fig, update, frames=n_points, interval=interval, blit=False, repeat=False)

    if args.save_animation:
        # Save the animation as an mp4 file
        ani.save('animation.mp4', writer='ffmpeg')

plt.show()




