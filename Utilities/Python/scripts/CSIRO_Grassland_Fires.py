
# Special script that plots flame fronts for CSIRO Grassland Fire Cases C064 and F19, fine.
# This script is not run by firebot. It was designed specifically to make a figure for a paper.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import struct
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../../exp/CSIRO_Grassland_Fires/'
outdir = '../../../../fds/Validation/CSIRO_Grassland_Fires/Current_Results/'
pltdir = '../../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CSIRO_Grassland_Fires/'

git_file = outdir + 'Case_C064_fine_cat_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


# Define classes to hold mesh data

class MeshData:
    def __init__(self):
        self.xs = None
        self.xf = None
        self.ys = None
        self.yf = None
        self.x = None
        self.y = None

class MeshGridData:
    def __init__(self):
        self.X = None
        self.Y = None

M = {}

# from Case_F19_fine input file
# &MESH IJK=80,80,80, XB=-20.0,0.0,0.0,20.0,0.0,20.0, MULT_ID='mesh' /
# &MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=11, J_LOWER=-6, J_UPPER=5, K_UPPER=0 /

XB = [-20.0, 0.0, 0.0, 20.0, 0.0, 20.0]
NX = 80
NY = 80
DX = 20
DY = 20
dx = DX / NX
dy = DY / NY
I_LOWER = 0
I_UPPER = 11
J_LOWER = -6
J_UPPER = 5

M3 = {}

NM3 = 0
for j in range(J_LOWER, J_UPPER + 1):
    for i in range(I_LOWER, I_UPPER + 1):
        NM3 = NM3 + 1
        M[NM3] = MeshData()
        M[NM3].xs = XB[0] + i * DX
        M[NM3].xf = XB[1] + i * DX
        M[NM3].ys = XB[2] + j * DY
        M[NM3].yf = XB[3] + j * DY
        M[NM3].x = np.arange(M[NM3].xs, M[NM3].xf + dx/2, dx)
        M[NM3].y = np.arange(M[NM3].ys, M[NM3].yf + dy/2, dy)
        M3[NM3] = MeshGridData()
        M3[NM3].X, M3[NM3].Y = np.meshgrid(M[NM3].x, M[NM3].y)

# from Case_C064_fine input file
# &MESH IJK=80,80,80, XB=-10.0,10.0,0.0,20.0,0.0,20.0, MULT_ID='mesh' /
# &MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=5, J_LOWER=-3, J_UPPER=2, K_UPPER=0 /

XB = [-10.0, 10.0, 0.0, 20.0, 0.0, 20.0]
NX = 80
NY = 80
DX = 20
DY = 20
dx = DX / NX
dy = DY / NY
I_LOWER = 0
I_UPPER = 5
J_LOWER = -3
J_UPPER = 2

M4 = {}

NM4 = 0
for j in range(J_LOWER, J_UPPER + 1):
    for i in range(I_LOWER, I_UPPER + 1):
        NM4 = NM4 + 1
        M[NM4] = MeshData()
        M[NM4].xs = XB[0] + i * DX
        M[NM4].xf = XB[1] + i * DX
        M[NM4].ys = XB[2] + j * DY
        M[NM4].yf = XB[3] + j * DY
        M[NM4].x = np.arange(M[NM4].xs, M[NM4].xf + dx/2, dx)
        M[NM4].y = np.arange(M[NM4].ys, M[NM4].yf + dy/2, dy)
        M4[NM4] = MeshGridData()
        M4[NM4].X, M4[NM4].Y = np.meshgrid(M[NM4].x, M[NM4].y)


def slread( fname, Tstart, Tend, *args, **kwargs ):

    """
    Reads FDS slice file
    (QQ,Time)=slread(fname,Tstart,Tend [,Nframes, gridskip, timeskip]);
      Tstart  is start time
      Tend    is end time
      Nframes is number of slice file frames (FDS default is 1000 between Tstart and Tend)
      gridskip is a skip rate for reading cells: =2 --> read every other cell, etc.
      timeskip is a skip rate for frames: =2 --> read every other frame, etc.
      QQ      contains the data
      Time    contains the time points
    """

    if len(args)==0:
        Nframes = 1000
    else:
        Nframes = args[0]
        
    gridskip = kwargs['gridskip'] if 'gridskip' in kwargs else 1
    timeskip = kwargs['timeskip'] if 'timeskip' in kwargs else 1
        
    f = open(fname,'rb')

    f.read(4)
    Str1 = f.read(30)       # quantity
    f.read(8)
    Str2 = f.read(30)       # short name
    f.read(8)
    Str3 = f.read(30)       # units
    f.read(8)
    Indx = struct.unpack('6i',f.read(24))  # index bounds i,j,k
    f.read(4)

    # allocate arrays for QQ and Time

    Isize = Indx[1]-Indx[0]+1
    Jsize = Indx[3]-Indx[2]+1
    Ksize = Indx[5]-Indx[4]+1
    if Isize==1:
       M = Jsize
       N = Ksize
    elif Jsize==1:
       M = Isize
       N = Ksize
    else:
       M = Isize
       N = Jsize

    Nframes = max(1,Nframes)
    QQ = np.zeros((M,N))
    Time = np.zeros(Nframes+1)
    
    ii = np.arange(0,M,gridskip)
    jj = np.arange(0,N,gridskip)
    tt = np.arange(0,Nframes+1,timeskip)
    Qskip  = np.zeros((len(ii),len(jj),len(tt)))

    st = 0

    while Time[st] < Tstart:

        f.read(4)
        Time[st] = np.fromfile(f, dtype=np.float32, count=1)[0]
        f.read(8)
        QQ = np.fromfile(f, dtype=np.float32, count=N*M).reshape(M, N)
        f.read(4)

    while Time[st] < Tend:

        st = st + 1
        f.read(4)
        Time[st] = np.fromfile(f, dtype=np.float32, count=1)[0]
        f.read(8)
        QQ = np.fromfile(f, dtype=np.float32, count=N*M).reshape(M, N)
        if st%timeskip==0:
            Qskip[:,:,int(st/timeskip)] = QQ[np.ix_(ii,jj)]
            Qskip[0,:,int(st/timeskip)] = 0.
            Qskip[:,0,int(st/timeskip)] = 0.
        f.read(4)
        if st>Nframes:
            break

    return(Qskip,Time[tt],int(st/timeskip))


# Case F19

fig = fdsplotlib.plot_to_fig(x_data=[-40,-40], y_data=[-40,-40],
                             x_min=-20, x_max=220, y_min=-120, y_max=120,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Windward Direction (m)', y_label='Crosswind Direction (m)')

T = np.array([56, 86, 138]) 
for k in range(len(T)):
    for n in range(1,NM3+1):
        nsf = 5
        if n >= 73 and n <= 84:
            nsf = 10
        Qdata, Time, nt = slread(outdir + 'Case_F19_fine_cat_' + str(n) + '_' + str(nsf) + '.sf', T[k], T[k] + 0.1)
        plt.contour(M3[n].X, M3[n].Y, Qdata[:, :,nt], levels=[20,100])

M = pd.read_csv(expdir + 'Case_F19_flame_position_data.csv', sep=',', header=0)

fl_56_x = M['56 x'].values
fl_56_y = M['56 y'].values
fl_86_x = M['86 x'].values
fl_86_y = M['86 y'].values
fl_138_x = M['138 x'].values
fl_138_y = M['138 y'].values

fdsplotlib.plot_to_fig(x_data=fl_56_x, y_data=fl_56_y, figure_handle=fig, marker_style='ko', data_label='56 s')
fdsplotlib.plot_to_fig(x_data=fl_86_x, y_data=fl_86_y, figure_handle=fig, marker_style='kv', data_label='86 s')
fdsplotlib.plot_to_fig(x_data=fl_138_x, y_data=fl_138_y, figure_handle=fig, marker_style='ks', data_label='138 s')

fdsplotlib.plot_to_fig(x_data=[0  ,200], y_data=[-100,-100], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[0  ,200], y_data=[ 100, 100], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[0  ,  0], y_data=[-100, 100], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[200,200], y_data=[-100, 100], figure_handle=fig, marker_style='k:')

plt.savefig(pltdir + 'Case_F19_flame_position.pdf', format='pdf')
plt.close()


# Case C064

fig = fdsplotlib.plot_to_fig(x_data=[-20,-20], y_data=[-20,-20],
                             x_min=-10, x_max=110, y_min=-60, y_max=60,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Windward Direction (m)', y_label='Crosswind Direction (m)')

T = np.array([27, 53, 100])
for k in range(len(T)):
    for n in range(1, NM4 + 1):
        nsf = 5
        if n >= 19 and n <= 24:
            nsf = 10
        Qdata, Time, nt = slread(outdir + 'Case_C064_fine_cat_' + str(n) + '_' + str(nsf) + '.sf', T[k], T[k] + 0.1)
        plt.contour(M4[n].X, M4[n].Y, Qdata[:, :,nt], levels=[20,100])

M = pd.read_csv(expdir + 'Case_C064_flame_position_data.csv', sep=',', header=0)

fl_27_x = M['27 x'].values
fl_27_y = M['27 y'].values
fl_53_x = M['53 x'].values
fl_53_y = M['53 y'].values
fl_100_x = M['100 x'].values
fl_100_y = M['100 y'].values

fdsplotlib.plot_to_fig(x_data=fl_27_x, y_data=fl_27_y, figure_handle=fig, marker_style='ko', data_label='56 s')
fdsplotlib.plot_to_fig(x_data=fl_53_x, y_data=fl_53_y, figure_handle=fig, marker_style='kv', data_label='86 s')
fdsplotlib.plot_to_fig(x_data=fl_100_x, y_data=fl_100_y, figure_handle=fig, marker_style='ks', data_label='138 s')

fdsplotlib.plot_to_fig(x_data=[0  ,100], y_data=[-50,-50], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[0  ,100], y_data=[ 50, 50], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[0  ,  0], y_data=[-50, 50], figure_handle=fig, marker_style='k:')
fdsplotlib.plot_to_fig(x_data=[100,100], y_data=[-50, 50], figure_handle=fig, marker_style='k:')

plt.savefig(pltdir + 'Case_C064_flame_position.pdf', format='pdf')
plt.close()

