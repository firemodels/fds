
# This script creates the Decaying Isotropic Turbulence plots in the FDS Verification Guide, chapter Turbulence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.patches import FancyArrowPatch
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../Verification/Turbulence/'
outdir = '../../Verification/Turbulence/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'csmag_32_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def plotspec_uvw(filename, marker_format, fig_handle):

    with open(filename, 'r') as f:
        header_line = f.readline().strip()
        data_lines = f.readlines()

    # Parse header to extract parameters
    s = [float(x) for x in header_line.split(',')]
    n = int(s[1])

    # Parse data into numpy array
    data = []
    for line in data_lines:
        if line.strip():  # Skip empty lines
            values = [float(x) for x in line.strip().split(',')]
            data.append(values)
    data = np.array(data)

    # Initialize 3D arrays
    u = np.zeros((n, n, n))
    v = np.zeros((n, n, n))
    w = np.zeros((n, n, n))
    tke = np.zeros((n, n, n))

    # Convert to 3D array
    p = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                u[i, j, k] = data[p, 0]
                v[i, j, k] = data[p, 1]
                w[i, j, k] = data[p, 2]
                tke[i, j, k] = 0.5 * (u[i, j, k]**2 + v[i, j, k]**2 + w[i, j, k]**2)
                p = p + 1

    # Perform FFT on velocity components
    u_hat = np.fft.fftn(u) / n**3
    v_hat = np.fft.fftn(v) / n**3
    w_hat = np.fft.fftn(w) / n**3

    # Compute spectral TKE
    tke_hat = np.zeros((n, n, n), dtype=complex)
    for k in range(n):
        for j in range(n):
            for i in range(n):
                tke_hat[i, j, k] = 0.5 * (u_hat[i, j, k] * np.conj(u_hat[i, j, k]) +
                                          v_hat[i, j, k] * np.conj(v_hat[i, j, k]) +
                                          w_hat[i, j, k] * np.conj(w_hat[i, j, k]))

    # spectrum
    L = 9 * 2 * np.pi / 100
    k0 = 2 * np.pi / L
    kmax = n // 2
    wn = k0 * np.arange(0, n + 1)  # wavenumber array
    vt = np.zeros(len(wn))

    # Bin spectral data into wavenumber bins
    for kx in range(1, n + 1):
        rkx = kx - 1
        if kx > kmax + 1:
            rkx = rkx - n

        for ky in range(1, n + 1):
            rky = ky - 1
            if ky > kmax + 1:
                rky = rky - n

            for kz in range(1, n + 1):
                rkz = kz - 1
                if kz > kmax + 1:
                    rkz = rkz - n

                rk = np.sqrt(rkx**2 + rky**2 + rkz**2)
                k = round(rk)

                vt[k] = vt[k] + tke_hat[kx - 1, ky - 1, kz - 1].real / k0

    # plot the energy spectrum
    wn_safe = np.where(wn[1:n] <= 0, 1e-9, wn[1:n])
    vt_safe = np.where(vt[1:n] <= 0, 1e-9, vt[1:n])
    fdsplotlib.plot_to_fig(x_data=wn_safe, y_data=vt_safe, figure_handle=fig_handle, marker_style=marker_format)


def energy_decay(chid, N):
    
    # Plot energy decay for Comte-Bellot and Corrsin data
    
    # chid = CHID from FDS input file
    # N = number of cells in 1D
    
    L = 0.56549  # \approx 9*(2*pi)/100;
    
    # Gather the Comte-Bellot/Corrsin data
    M = np.loadtxt(expdir + 'cbcdata.txt')
    k = M[:, 0] * 1e2
    E1 = M[:, 1] / 1e6
    E2 = M[:, 2] / 1e6
    E3 = M[:, 3] / 1e6
    
    # Apply transfer function
    k0 = 2 * np.pi / L
    kc = 0.5 * N * k0
    delta = np.pi / kc
    # G = sin(1/2*k.*delta)./(1/2*k.*delta);  # box filter
    G = np.ones(len(k))
    for j in range(len(k)):
        if k[j] > kc:
            tmp = 3 * kc / k[j] - 2
            if tmp > 0:
                G[j] = np.sqrt(tmp)  # RJM filter
            else:
                G[j] = 0
    
    E1 = G * G * E1
    E2 = G * G * E2
    E3 = G * G * E3
    
    # Now integrate
    E1_bar = 0
    E2_bar = 0
    E3_bar = 0
    for j in range(len(k) - 1):
        dk = k[j + 1] - k[j]
        # E1_bar = E1_bar + sqrt(E1[j]*E1[j+1])*dk
        # E2_bar = E2_bar + sqrt(E2[j]*E2[j+1])*dk
        # E3_bar = E3_bar + sqrt(E3[j]*E3[j+1])*dk
        E1_bar = E1_bar + 0.5 * (E1[j] + E1[j + 1]) * dk
        E2_bar = E2_bar + 0.5 * (E2[j] + E2[j + 1]) * dk
        E3_bar = E3_bar + 0.5 * (E3[j] + E3[j + 1]) * dk
    
    # Create figure
    fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=0.7, y_min=0, y_max=0.07,
                             revision_label=version_string,
                             x_label='Time (s)',
                             y_label='Kinetic Energy (m$^2$/s$^2$)')

    if chid == 'csmag_32':
        M = pd.read_csv(outdir + 'mu0_32_devc.csv', skiprows=2, header=None).values
        t = M[:, 0]
        ke = M[:, 1]
        fdsplotlib.plot_to_fig(x_data=t, y_data=ke, figure_handle=fig, marker_style='k-.', data_label='FDS zero visc')
        
        M = pd.read_csv(outdir + 'csmag0_32_devc.csv', skiprows=2, header=None).values
        t = M[:, 0]
        ke = M[:, 1]
        fdsplotlib.plot_to_fig(x_data=t, y_data=ke, figure_handle=fig, marker_style='k--', data_label='FDS mol visc')
    
    M = pd.read_csv(outdir + chid + '_devc.csv', skiprows=2, header=None).values
    t = M[:, 0]
    ke = M[:, 1]
    if chid=='csmag_32': 
        label='FDS eddy visc'
    else:
        label='FDS'
    fdsplotlib.plot_to_fig(x_data=t, y_data=ke, figure_handle=fig, marker_style='k-', data_label=label)
    
    # Plot filtered CBC data points
    fdsplotlib.plot_to_fig(x_data=0.00, y_data=E1_bar, figure_handle=fig, marker_style='ko', data_label='Filtered CBC data')
    fdsplotlib.plot_to_fig(x_data=0.28, y_data=E2_bar, figure_handle=fig, marker_style='ro')
    fdsplotlib.plot_to_fig(x_data=0.66, y_data=E3_bar, figure_handle=fig, marker_style='bo')
    
    plt.savefig(pltdir + chid + '_decay.pdf', format='pdf')
    plt.close()


def plotspec(chid, N):

    # Read and plot energy spectra for Comte-Bellot and Corrsin data

    # chid: CHID from FDS input file
    # N  Number of cells in 1D

    fig = fdsplotlib.plot_to_fig(x_data=[1e-10,1e-10], y_data=[1e-10,1e-10],
                                 x_min=10, x_max=2000, y_min=1e-7, y_max=1e-3,
                                 revision_label=version_string,
                                 plot_type='loglog',
                                 x_label='$k$ (1/m)',
                                 y_label='$E(k)$ (m$^3$/s$^2$)')

    L = 9*2*np.pi/100  # box length (m)
    k0 = 2*np.pi/L
    kc = 1/2*N*k0

    uvw_file1 = f"{outdir}{chid}_uvw_t0_m1.csv"
    uvw_file2 = f"{outdir}{chid}_uvw_t1_m1.csv"
    uvw_file3 = f"{outdir}{chid}_uvw_t2_m1.csv"

    # Plot the FDS data
    plotspec_uvw(uvw_file1, 'ko-', fig)
    plotspec_uvw(uvw_file2, 'ro-', fig)
    plotspec_uvw(uvw_file3, 'bo-', fig)

    # Gather the Comte-Bellot/Corrsin data
    CBC = np.loadtxt(f'{outdir}cbcdata.txt')
    k  = CBC[1:, 0] * 1e2
    E1 = CBC[1:, 1] / 1e6
    E2 = CBC[1:, 2] / 1e6
    E3 = CBC[1:, 3] / 1e6

    # Plot the CBC data
    fdsplotlib.plot_to_fig(x_data=k, y_data=E1, figure_handle=fig, marker_style='k-')
    fdsplotlib.plot_to_fig(x_data=k, y_data=E2, figure_handle=fig, marker_style='k-')
    fdsplotlib.plot_to_fig(x_data=k, y_data=E3, figure_handle=fig, marker_style='k-')

    # Add lines to indicate cutoff wavenumbers
    fdsplotlib.plot_to_fig(x_data=[k0*N/2, k0*N/2], y_data=[1e-10, 1e-2], figure_handle=fig, marker_style='k--')  # LES Nyquist limit

    # Add annotation arrow
    ax = plt.gca()
    ax.annotate('Time',
                xy=(0.44, 0.55), xycoords='axes fraction',
                xytext=(0.56, 0.85), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='black'))

    # # Apply transfer function (optional)
    # delta = np.pi / kc
    # G = np.zeros(len(k))
    # for j in range(len(k)):
    #     if k[j] <= kc:
    #         # G[j] = np.sin(1/2*k[j]*delta) / (1/2*k[j]*delta)  # box filter
    #         # G[j] = np.exp(-k[j]**2*delta**2/24)                # Gaussian filter
    #         G[j] = 1                                             # spectral cutoff
    #     elif k[j] > kc and k[j] <= np.sqrt(2)*kc:
    #         G[j] = np.sqrt(3*kc/k[j] - 2)                        # RJM implicit filter
    #     elif k[j] > np.sqrt(2)*kc:
    #         G[j] = 0
    #     E1[j] = G[j]*G[j]*E1[j]
    #     E2[j] = G[j]*G[j]*E2[j]
    #     E3[j] = G[j]*G[j]*E3[j]
    #
    # # Plot the filtered CBC data
    # plt.loglog(k, E1, 'ko', linewidth=2)
    # plt.loglog(k, E2, 'ro', linewidth=2)
    # plt.loglog(k, E3, 'bo', linewidth=2)

    plt.savefig(f'../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/{chid}_spectra.pdf', format='pdf')
    plt.close()

energy_decay('csmag_32',32)
energy_decay('csmag_64',64)
energy_decay('dsmag_32',32)
energy_decay('dsmag_64',64)
energy_decay('deardorff_32',32)
energy_decay('deardorff_64',64)
energy_decay('vreman_32',32)
energy_decay('vreman_64',64)
energy_decay('wale_32',32)
energy_decay('wale_64',64)

plotspec('csmag_32',32)
plotspec('csmag_64',64)
plotspec('dsmag_32',32)
plotspec('dsmag_64',64)
plotspec('deardorff_32',32)
plotspec('deardorff_64',64)
plotspec('vreman_32',32)
plotspec('vreman_64',64)
plotspec('wale_32',32)
plotspec('wale_64',64)

