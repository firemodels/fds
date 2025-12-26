
# Particle/Gas momentum transfer

import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Sprinklers_and_Sprays/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'fluid_part_mom_x_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def read_prt5(filename, precision, *varargin):
    """
    Read FDS part file (*.prt5)

    Parameters:
    -----------
    filename : str
        Name of the .prt5 file
    precision : str
        'real*4' or 'real*8' depending on EB_PART_FILE setting
    *varargin : optional arguments
        Pass 'evac', 'EVAC', or 'Evac' for FDS+Evac files

    Returns:
    --------
    STIME : numpy.ndarray
        Simulation times
    XP : numpy.ndarray
        X positions of particles
    YP : numpy.ndarray
        Y positions of particles
    ZP : numpy.ndarray
        Z positions of particles
    QP : numpy.ndarray
        Particle quantities
    AP : numpy.ndarray (optional, only for evac mode with 6 outputs)
        Additional particle data for evacuation mode
    """

    # Determine number of outputs requested by caller
    # In Python, we'll return a tuple and let caller unpack as needed
    # For now, we'll track if AP should be returned
    nout = 6 if len(varargin) > 0 else 5

    evac = False
    if len(varargin) > 0:
        if varargin[0] in ['evac', 'EVAC', 'Evac']:
            evac = True

    fid = open(filename, 'rb')

    # Determine precision formats for struct
    if precision == 'real*4':
        real_format = '<f'  # little-endian float (4 bytes)
        real_dtype = np.float32
        real_size = 4
    elif precision == 'real*8':
        real_format = '<d'  # little-endian double (8 bytes)
        real_dtype = np.float64
        real_size = 8
    else:
        raise ValueError("precision must be 'real*4' or 'real*8'")

    int_format = '<i'  # little-endian integer (4 bytes)
    int_size = 4

    # The DUMMY lines are 4 byte placeholders that apparently fortran puts at the
    # beginning and end of all lines.  I only knew of this thanks to Glenn.

    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
    ONE_INTEGER = struct.unpack(int_format, fid.read(int_size))[0]
    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
    INT_VERSION = struct.unpack(int_format, fid.read(int_size))[0]
    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
    N_PART = struct.unpack(int_format, fid.read(int_size))[0]
    DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

    N_QUANTITIES = []
    SMOKEVIEW_LABEL = []
    UNITS = []

    for NPC in range(N_PART):
        DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
        PC = struct.unpack('<2i', fid.read(2 * int_size))
        N_QUANTITIES.append(PC[0])
        DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

        smokeview_labels_npc = []
        units_npc = []

        for NQ in range(N_QUANTITIES[NPC]):
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
            # Read 30 bytes for label
            label_bytes = fid.read(30)
            # Decode and strip null characters and whitespace
            label = label_bytes.decode('ascii', errors='ignore').rstrip('\x00').rstrip()
            smokeview_labels_npc.append(label)
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
            # Read 30 bytes for units
            units_bytes = fid.read(30)
            # Decode and strip null characters and whitespace
            unit = units_bytes.decode('ascii', errors='ignore').rstrip('\x00').rstrip()
            units_npc.append(unit)
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

        SMOKEVIEW_LABEL.append(smokeview_labels_npc)
        UNITS.append(units_npc)

    # Initialize lists to store data (will convert to arrays later)
    STIME_list = []
    XP_list = []
    YP_list = []
    ZP_list = []
    QP_list = []
    AP_list = []

    n = 0
    while True:
        # Try to read DUMMY
        dummy_bytes = fid.read(int_size)
        if len(dummy_bytes) < int_size:
            break
        DUMMY = struct.unpack(int_format, dummy_bytes)[0]

        # Try to read stime
        stime_bytes = fid.read(real_size)
        if len(stime_bytes) < real_size:
            break
        stime_tmp = struct.unpack(real_format, stime_bytes)[0]

        DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

        if len(stime_bytes) == 0:
            break
        else:
            STIME_list.append(stime_tmp)

        # Initialize temporary storage for this time step
        XP_time = []
        YP_time = []
        ZP_time = []
        QP_time = []
        AP_time = []

        for NPC in range(N_PART):
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
            NPLIM = struct.unpack(int_format, fid.read(int_size))[0]
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

            # Read particle positions
            xp_bytes = fid.read(NPLIM * real_size)
            xp = np.frombuffer(xp_bytes, dtype=real_dtype)

            yp_bytes = fid.read(NPLIM * real_size)
            yp = np.frombuffer(yp_bytes, dtype=real_dtype)

            zp_bytes = fid.read(NPLIM * real_size)
            zp = np.frombuffer(zp_bytes, dtype=real_dtype)

            if evac:
                ap1_bytes = fid.read(NPLIM * real_size)
                ap1 = np.frombuffer(ap1_bytes, dtype=real_dtype)

                ap2_bytes = fid.read(NPLIM * real_size)
                ap2 = np.frombuffer(ap2_bytes, dtype=real_dtype)

                ap3_bytes = fid.read(NPLIM * real_size)
                ap3 = np.frombuffer(ap3_bytes, dtype=real_dtype)

                ap4_bytes = fid.read(NPLIM * real_size)
                ap4 = np.frombuffer(ap4_bytes, dtype=real_dtype)

            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

            # Store particle positions for this class
            XP_time.append(xp)
            YP_time.append(yp)
            ZP_time.append(zp)

            if evac and nout == 6:
                # Stack the four additional parameters
                ap_combined = np.stack([ap1, ap2, ap3, ap4], axis=-1)
                AP_time.append(ap_combined)

            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
            TA_bytes = fid.read(NPLIM * int_size)
            TA = np.frombuffer(TA_bytes, dtype=np.int32)
            DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

            if N_QUANTITIES[NPC] > 0:
                DUMMY = struct.unpack(int_format, fid.read(int_size))[0]

                qp_class = []
                for NQ in range(N_QUANTITIES[NPC]):
                    qp_bytes = fid.read(NPLIM * real_size)
                    qp = np.frombuffer(qp_bytes, dtype=real_dtype)
                    qp_class.append(qp)

                QP_time.append(qp_class)
                DUMMY = struct.unpack(int_format, fid.read(int_size))[0]
            else:
                QP_time.append([])

        XP_list.append(XP_time)
        YP_list.append(YP_time)
        ZP_list.append(ZP_time)
        QP_list.append(QP_time)
        if evac and nout == 6:
            AP_list.append(AP_time)

        n += 1

    fid.close()

    # Convert lists to numpy arrays with proper dimensions
    STIME = np.array(STIME_list)

    # Determine maximum number of particles across all time steps and classes
    max_particles = 0
    for time_step in XP_list:
        for particle_class in time_step:
            max_particles = max(max_particles, len(particle_class))

    # Initialize arrays with NaN for missing data
    n_time = len(STIME_list)
    XP = np.full((n_time, max_particles, N_PART), np.nan)
    YP = np.full((n_time, max_particles, N_PART), np.nan)
    ZP = np.full((n_time, max_particles, N_PART), np.nan)

    # Fill position arrays
    for t in range(n_time):
        for npc in range(N_PART):
            n_particles = len(XP_list[t][npc])
            XP[t, :n_particles, npc] = XP_list[t][npc]
            YP[t, :n_particles, npc] = YP_list[t][npc]
            ZP[t, :n_particles, npc] = ZP_list[t][npc]

    # Handle QP array (4D)
    max_quantities = max(N_QUANTITIES) if N_QUANTITIES else 0
    QP = np.full((n_time, max_particles, N_PART, max_quantities), np.nan)

    for t in range(n_time):
        for npc in range(N_PART):
            if len(QP_list[t][npc]) > 0:
                n_particles = len(QP_list[t][npc][0])
                for nq in range(len(QP_list[t][npc])):
                    QP[t, :n_particles, npc, nq] = QP_list[t][npc][nq]

    # Handle AP array for evac mode
    if evac and nout == 6:
        AP = np.full((n_time, max_particles, N_PART, 4), np.nan)
        for t in range(n_time):
            for npc in range(N_PART):
                n_particles = AP_list[t][npc].shape[0]
                AP[t, :n_particles, npc, :] = AP_list[t][npc]

        return STIME, XP, YP, ZP, QP, AP
    else:
        return STIME, XP, YP, ZP, QP


# Read CSV files (skip 2 header rows)
M = pd.read_csv(outdir + 'fluid_part_mom_x_devc.csv', skiprows=1)
tx = M.iloc[:, 0].values
U = M.iloc[:, 1].values
MX = M.iloc[:, 2].values

M = pd.read_csv(outdir + 'fluid_part_mom_y_devc.csv', skiprows=1)
ty = M.iloc[:, 0].values
V = M.iloc[:, 1].values
MY = M.iloc[:, 2].values

M = pd.read_csv(outdir + 'fluid_part_mom_z_devc.csv', skiprows=1)
tz = M.iloc[:, 0].values
W = M.iloc[:, 1].values
MZ = M.iloc[:, 2].values

# Create range for plotting (every 10th point)
range_idx = np.arange(0, min([len(tx), len(ty), len(tz)]), 10)

# Create first figure
fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=1, y_min=0, y_max=150,
                             revision_label=version_string,
                             figure_size=(plot_style['Paper_Width']+0.4,plot_style['Paper_Height']),
                             legend_location='outside',
                             x_label='Time (s)',
                             y_label='Momentum (kg m/s)')

fdsplotlib.plot_to_fig(x_data=tx[range_idx], y_data=MX[range_idx]*U[range_idx], figure_handle=fig, marker_style='bo', marker_fill_color='none', data_label='FDS fluid $u$')
fdsplotlib.plot_to_fig(x_data=ty[range_idx], y_data=MY[range_idx]*V[range_idx], figure_handle=fig, marker_style='bv', marker_fill_color='none', data_label='FDS fluid $v$')
fdsplotlib.plot_to_fig(x_data=tz[range_idx], y_data=MZ[range_idx]*W[range_idx], figure_handle=fig, marker_style='b+', marker_fill_color='none', data_label='FDS fluid $w$')

# Particle parameters
n = 1000                  # number of particles
mpv = 10.                 # kg/m^3, mass_per_volume (from FDS input file)
v_xb = 1**3               # volume of XB region on init line in FDS input file
rho_p = 1000              # density of water, kg/m^3
d_p = 1000e-6             # diameter, m
v_p = 4/3 * np.pi * (d_p/2)**3   # volume of a single droplet, m^3
m_p = rho_p * v_p         # mass of single droplet, kg
pwt = mpv * v_xb / (n * m_p)     # particle weight factor

# Read particle data from PRT5 files
STIME_X, XP_X, YP_X, ZP_X, QP_X = read_prt5(outdir + 'fluid_part_mom_x_1.prt5', 'real*4')
STIME_Y, XP_Y, YP_Y, ZP_Y, QP_Y = read_prt5(outdir + 'fluid_part_mom_y_1.prt5', 'real*4')
STIME_Z, XP_Z, YP_Z, ZP_Z, QP_Z = read_prt5(outdir + 'fluid_part_mom_z_1.prt5', 'real*4')

# Calculate particle momentum
P_X = np.zeros(len(STIME_X))
for i in range(len(STIME_X)):
    for j in range(n):
        P_X[i] = P_X[i] + pwt * m_p * QP_X[i, j, 0, 0]  # momentum of particle at time STIME(i)

P_Y = np.zeros(len(STIME_Y))
for i in range(len(STIME_Y)):
    for j in range(n):
        P_Y[i] = P_Y[i] + pwt * m_p * QP_Y[i, j, 0, 0]

P_Z = np.zeros(len(STIME_Z))
for i in range(len(STIME_Z)):
    for j in range(n):
        P_Z[i] = P_Z[i] + pwt * m_p * QP_Z[i, j, 0, 0]

# Create range for particle plotting
range_p = np.arange(0, min([len(STIME_X), len(STIME_Y), len(STIME_Z)]), 10)

# Plot particle momentum
fdsplotlib.plot_to_fig(x_data=STIME_X[range_p], y_data=P_X[range_p], figure_handle=fig, marker_style='ro', marker_fill_color='none', data_label=r'FDS particle $u_\mathrm{p}$')
fdsplotlib.plot_to_fig(x_data=STIME_Y[range_p], y_data=P_Y[range_p], figure_handle=fig, marker_style='rv', marker_fill_color='none', data_label=r'FDS particle $v_\mathrm{p}$')
fdsplotlib.plot_to_fig(x_data=STIME_Z[range_p], y_data=P_Z[range_p], figure_handle=fig, marker_style='r+', marker_fill_color='none', data_label=r'FDS particle $w_\mathrm{p}$')

# Calculate total momentum
P_total_X = P_X + (MX * U)  # total momentum
P_total_Y = P_Y + (MY * V)
P_total_Z = P_Z + (MZ * W)

# Plot total momentum
fdsplotlib.plot_to_fig(x_data=STIME_X[range_p], y_data=P_total_X[range_p], figure_handle=fig, marker_style='go', marker_fill_color='none', data_label=r'FDS total $U$')
fdsplotlib.plot_to_fig(x_data=STIME_Y[range_p], y_data=P_total_Y[range_p], figure_handle=fig, marker_style='gv', marker_fill_color='none', data_label=r'FDS total $V$')
fdsplotlib.plot_to_fig(x_data=STIME_Z[range_p], y_data=P_total_Z[range_p], figure_handle=fig, marker_style='g+', marker_fill_color='none', data_label=r'FDS total $W$')

# analytical solution
rho = 1.1992661     # fluid density
Cd = 1              # drag coefficient
A_p = np.pi * (d_p/2)**2  # particle area

u_p = QP_Z[0, 0, 0, 0]  # initial velocity
U_p = U[0]

M_p = MX[0] / (n * pwt)
alpha = M_p / m_p

# Initialize solution arrays
u_soln = []
U_soln = []
t_soln = []

for i in range(len(tx) - 1):
    u_soln.append(u_p)
    U_soln.append(U_p)
    t_soln.append(tx[i])

    dt = tx[i + 1] - tx[i]

    u0 = u_p
    U0 = U_p

    beta = 0.5 * rho * Cd * A_p * (1/m_p + 1/M_p) * abs(u0 - U0)
    u_p = u0 / (1 + beta * dt) + (u0 + alpha * U0) / (1 + alpha) * (beta * dt) / (1 + beta * dt)
    U_p = U0 + n * pwt * m_p / MX[0] * (u0 - u_p)

# Convert lists to arrays
u_soln = np.array(u_soln)
U_soln = np.array(U_soln)
t_soln = np.array(t_soln)

# Plot analytical solutions
fdsplotlib.plot_to_fig(x_data=t_soln, y_data=MX[0]*U_soln,     figure_handle=fig, marker_style='b-', data_label='Analytical fluid')
fdsplotlib.plot_to_fig(x_data=t_soln, y_data=n*pwt*m_p*u_soln, figure_handle=fig, marker_style='r-', data_label='Analytical particle')

plt.savefig(pltdir + 'fluid_part_momentum.pdf', format='pdf')
plt.close()


# plot velocities

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=1, y_min=0, y_max=10,
                             revision_label=version_string,
                             figure_size=(plot_style['Paper_Width']+0.4,plot_style['Paper_Height']),
                             legend_location='outside',
                             x_label='Time (s)',
                             y_label='Velocity (m/s)')

# Plot fluid velocities
fdsplotlib.plot_to_fig(x_data=tx[range_idx], y_data=U[range_idx], figure_handle=fig, marker_style='bo', marker_fill_color='none', data_label='FDS fluid $u$')
fdsplotlib.plot_to_fig(x_data=ty[range_idx], y_data=V[range_idx], figure_handle=fig, marker_style='bv', marker_fill_color='none', data_label='FDS fluid $v$')
fdsplotlib.plot_to_fig(x_data=tz[range_idx], y_data=W[range_idx], figure_handle=fig, marker_style='b+', marker_fill_color='none', data_label='FDS fluid $w$')

# Calculate particle velocities
U_p = P_X / (n * pwt * m_p)
V_p = P_Y / (n * pwt * m_p)
W_p = P_Z / (n * pwt * m_p)

# Plot particle velocities
fdsplotlib.plot_to_fig(x_data=tx[range_idx], y_data=U_p[range_idx], figure_handle=fig, marker_style='ro', marker_fill_color='none', data_label=r'FDS particle $u_\mathrm{p}$')
fdsplotlib.plot_to_fig(x_data=ty[range_idx], y_data=V_p[range_idx], figure_handle=fig, marker_style='rv', marker_fill_color='none', data_label=r'FDS particle $v_\mathrm{p}$')
fdsplotlib.plot_to_fig(x_data=tz[range_idx], y_data=W_p[range_idx], figure_handle=fig, marker_style='r+', marker_fill_color='none', data_label=r'FDS particle $w_\mathrm{p}$')

# Calculate equilibrium velocities
U_eq = P_total_X / (MX[0] + n * pwt * m_p)
V_eq = P_total_Y / (MY[0] + n * pwt * m_p)
W_eq = P_total_Z / (MZ[0] + n * pwt * m_p)

# Plot equilibrium velocities
fdsplotlib.plot_to_fig(x_data=tx, y_data=U_eq, figure_handle=fig, marker_style='g--', data_label='Equilibrium Velocity')
fdsplotlib.plot_to_fig(x_data=ty, y_data=V_eq, figure_handle=fig, marker_style='g--')
fdsplotlib.plot_to_fig(x_data=tz, y_data=W_eq, figure_handle=fig, marker_style='g--')

# Save figure
plt.savefig(pltdir + 'fluid_part_velocity.pdf', format='pdf')
plt.close()


# Particle Drag Profile (formerly part_drag_profile.m)

r_p = 0.001               # radius, m
l_p = 0.02                # length, m
v_p = np.pi*(r_p)**2*l_p  # volume of a single particle, m^3
shape_factor = 0.25                      # assumes random orientation of cylinders
a_p = shape_factor*l_p*(2*np.pi*r_p)     # projected area, m^2

z   = np.linspace(0, 10, 20)
u_z = z
c_d = 2.8        # from FDS input file (specified)
rho_g = 1.195    # from FDS out file
f_x = c_d * a_p * 0.5*rho_g*(u_z**2)  # drag experienced by a single particle

fig = fdsplotlib.plot_to_fig(x_data=z, y_data=-f_x, line_style='k-', data_label='exact',
                             x_min=0, x_max=10, y_min=-0.01, y_max=0,
                             revision_label=version_string,
                             x_label='Position (m)', y_label='Drag Force (N)',
                             legend_location='lower left')

ddir = '../../Verification/WUI/'
chid = ['part_drag_prof_ux', 'part_drag_prof_uy', 'part_drag_prof_uz',
        'part_drag_prof_vx', 'part_drag_prof_vy', 'part_drag_prof_vz',
        'part_drag_prof_wx', 'part_drag_prof_wy', 'part_drag_prof_wz']
j = [1, 2, 3, 1, 2, 3, 1, 2, 3]  # coordinate direction (x=1, y=2, z=3)

for i in range(len(chid)): 

    STIME, XP, YP, ZP, QP = read_prt5(ddir + chid[i] + '_1.prt5', 'real*4')

    N = 100

    if j[i] == 1:
        npts = XP.shape[1]
        idx = np.linspace(0, npts - 1, N).astype(int)
        fds_key_label = 'FDS part $x$' if i==0 else None
        fdsplotlib.plot_to_fig(x_data=XP[-1, idx], y_data=QP[-1, idx, 0, 0]/QP[-1, idx, 0, 1], figure_handle=fig, marker_style='bo', data_label=fds_key_label)
        v = np.abs(c_d * a_p * 0.5*rho_g*(XP[-1, idx]**2) - QP[-1, idx, 0, 0]/QP[-1, idx, 0, 1])
    elif j[i] == 2:
        npts = YP.shape[1]
        idy = np.linspace(0, npts - 1, N).astype(int)
        fds_key_label = 'FDS part $y$' if i==1 else None
        fdsplotlib.plot_to_fig(x_data=YP[-1, idy], y_data=QP[-1, idy, 0, 0]/QP[-1, idy, 0, 1], figure_handle=fig, marker_style='ro', data_label=fds_key_label)
        v = np.abs(c_d * a_p * 0.5*rho_g*(YP[-1, idy]**2) - QP[-1, idy, 0, 0]/QP[-1, idy, 0, 1])
    elif j[i] == 3:
        npts = ZP.shape[1]
        idz = np.linspace(0, npts - 1, N).astype(int)
        fds_key_label = 'FDS part $z$' if i==2 else None
        fdsplotlib.plot_to_fig(x_data=ZP[-1, idz], y_data=QP[-1, idz, 0, 0]/QP[-1, idz, 0, 1], figure_handle=fig, marker_style='go', data_label=fds_key_label)
        v = np.abs(c_d * a_p * 0.5*rho_g*(ZP[-1, idz]**2) - QP[-1, idz, 0, 0]/QP[-1, idz, 0, 1])

    err = np.linalg.norm(v)/len(v)
    #if err > 1e-4:
    #    print('Error: Case ' + ddir + chid[i] + ' error = ' + str(err))

plt.savefig(pltdir + 'part_drag_profile.pdf', format='pdf')
plt.close()


# Particle Temperature Profile

outdir = '../../Verification/WUI/'

git_file = outdir + 'part_temp_prof_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

ddir = '../../Verification/WUI/'
chid = 'part_temp_prof'

STIME, XP, YP, ZP, QP = read_prt5(ddir + chid + '_1.prt5', 'real*4')
N = 200
npts = ZP.shape[1]
idz = np.linspace(0, npts - 1, N).astype(int)
z_fds = ZP[-1, idz].flatten()
T_g_fds = QP[-1, idz, 0, 0].flatten()
fig = fdsplotlib.plot_to_fig(x_data=z_fds, y_data=T_g_fds, marker_style='bo',
                            x_min=0, x_max=10, y_min=0, y_max=1000,
                            data_label='FDS',revision_label=version_string,
                            x_label='Position (m)', y_label='Temperature (°C)')

z   = np.linspace(0, 10, 20)
T_g = 20+7.8*z**2
fdsplotlib.plot_to_fig(x_data=z, y_data=T_g, line_style='k-', data_label='exact',figure_handle=fig)

err = np.mean(np.abs(1-T_g_fds/(20+7.8*z_fds**2)))
if err > 5e-2:
   print('Error: Case ' + ddir + chid + ' error = ' + str(err))

plt.savefig(pltdir + 'part_temp_profile.pdf', format='pdf')
plt.close()

# Particle Species Profile

outdir = '../../Verification/WUI/'

git_file = outdir + 'part_spec_prof_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

ddir = '../../Verification/WUI/'
chid = 'part_spec_prof'

STIME, XP, YP, ZP, QP = read_prt5(ddir + chid + '_1.prt5', 'real*4')
N = 200
npts = ZP.shape[1]
idz = np.linspace(0, npts - 1, N).astype(int)
z_fds = ZP[-1, idz].flatten()
Y_O2_fds = QP[-1, idz, 0, 0].flatten()/QP[-1, idz, 0, 1].flatten()
fig = fdsplotlib.plot_to_fig(x_data=z_fds, y_data=Y_O2_fds, marker_style='bo',
                            x_min=0, x_max=10, y_min=0, y_max=0.23,
                            data_label='FDS',revision_label=version_string,
                            x_label='Position (m)', y_label='Oxygen Mass Fraction (-)')

z   = np.linspace(0, 10, 20)
Y_O2 = 0.0023*z**2
fdsplotlib.plot_to_fig(x_data=z, y_data=Y_O2, line_style='k-', data_label='exact',figure_handle=fig)

# Clip small values to not overly bias error
z_fds = z_fds[Y_O2_fds>0.001]
Y_O2_fds = Y_O2_fds[Y_O2_fds>0.001]

err = np.mean(np.abs(1-Y_O2_fds/(0.0023*z_fds**2)))
if err > 1e-1:
   print('Error: Case ' + ddir + chid + ' error = ' + str(err))

plt.savefig(pltdir + 'part_spec_profile.pdf', format='pdf')
plt.close()

