
# Solid body rotation (soborot) flow field, scalar transport

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Scalar_Analytical_Solution/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'soborot_charm_cos_wave_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

chid = ['soborot_charm_cos_wave_16',
        'soborot_charm_cos_wave_32',
        'soborot_charm_cos_wave_64',
        'soborot_charm_square_wave_128',
        'soborot_charm_square_wave_16',
        'soborot_charm_square_wave_32',
        'soborot_charm_square_wave_64',
        'soborot_godunov_square_wave_128',
        'soborot_godunov_square_wave_16',
        'soborot_godunov_square_wave_32',
        'soborot_godunov_square_wave_64',
        'soborot_superbee_cos_wave_128',
        'soborot_superbee_cos_wave_16',
        'soborot_superbee_cos_wave_32',
        'soborot_superbee_cos_wave_64',
        'soborot_superbee_square_wave_128',
        'soborot_superbee_square_wave_128_1mesh',
        'soborot_superbee_square_wave_16',
        'soborot_superbee_square_wave_32',
        'soborot_superbee_square_wave_64']

# positions of DEVC along diagonal D
L = 1.0
D = L * np.sqrt(2.0)
r_exact = np.linspace(0.0, D, 1000)
Y_exact = np.zeros_like(r_exact)
i_range = np.where((r_exact > 0.25) & (r_exact < 0.75))[0]
Y_exact[i_range] = 1.0
r_16 = np.arange((D/32.0), D, (D/16.0))
r_32 = np.arange((D/64.0), D, (D/32.0))
r_64 = np.arange((D/128.0), D, (D/64.0))
r_128 = np.arange((D/256.0), D, (D/128.0))

def read_tracers(csv_path, n_tracers):
    # csv has a header row followed by column headers in the 2nd row in Matlab importdata(...,2)
    # here assume the first row is a header and the next row are column names -> skiprows=0 to read header
    # Many CSVs from FDS include an initial line; using skiprows=1 like in previous conversion
    df = pd.read_csv(csv_path, skiprows=1)
    # strip whitespace from column names
    df.columns = df.columns.str.strip()
    col_start = df.columns.get_loc('Y_TRACER-1')
    col_end = df.columns.get_loc(f'Y_TRACER-{n_tracers}')
    return df.iloc[-1, col_start:col_end+1].values

# SQUARE WAVE: CHARM

Y_charm_16 = read_tracers(os.path.join(outdir, 'soborot_charm_square_wave_16_devc.csv'), 16)
Y_charm_32 = read_tracers(os.path.join(outdir, 'soborot_charm_square_wave_32_devc.csv'), 32)
Y_charm_64 = read_tracers(os.path.join(outdir, 'soborot_charm_square_wave_64_devc.csv'), 64)
Y_charm_128= read_tracers(os.path.join(outdir, 'soborot_charm_square_wave_128_devc.csv'),128)

fig = fdsplotlib.plot_to_fig(x_data=r_exact, y_data=Y_exact, data_label='Exact', line_style='k-',
                  x_label='Radial Position (m)', y_label='Scalar Mass Fraction',
                  plot_title='CHARM',
                  x_min=0, x_max=1, y_min=0, y_max=1.2)
fdsplotlib.plot_to_fig(x_data=r_16, y_data=Y_charm_16, figure_handle=fig, data_label='$n=16$', line_style='bo-')
fdsplotlib.plot_to_fig(x_data=r_32, y_data=Y_charm_32, figure_handle=fig, data_label='$n=32$', line_style='m*-')
fdsplotlib.plot_to_fig(x_data=r_64, y_data=Y_charm_64, figure_handle=fig, data_label='$n=64$', line_style='r^-')
fdsplotlib.plot_to_fig(x_data=r_128, y_data=Y_charm_128, figure_handle=fig, data_label='$n=128$', line_style='gsq-')

plt.savefig(os.path.join(pltdir, 'soborot_charm_square_wave.pdf'), format='pdf')
plt.close()

# SQUARE WAVE: SUPERBEE 

Y_superbee_16 = read_tracers(os.path.join(outdir, 'soborot_superbee_square_wave_16_devc.csv'), 16)
Y_superbee_32 = read_tracers(os.path.join(outdir, 'soborot_superbee_square_wave_32_devc.csv'), 32)
Y_superbee_64 = read_tracers(os.path.join(outdir, 'soborot_superbee_square_wave_64_devc.csv'), 64)
Y_superbee_128 = read_tracers(os.path.join(outdir, 'soborot_superbee_square_wave_128_devc.csv'), 128)
Y_superbee_128_1mesh = read_tracers(os.path.join(outdir, 'soborot_superbee_square_wave_128_1mesh_devc.csv'), 128)

fig = fdsplotlib.plot_to_fig(x_data=r_exact, y_data=Y_exact, data_label='Exact', line_style='k-',
                  x_label='Radial Position (m)', y_label='Scalar Mass Fraction',
                  plot_title='Superbee',
                  x_min=0, x_max=1, y_min=0, y_max=1.2)
fdsplotlib.plot_to_fig(x_data=r_16, y_data=Y_superbee_16, figure_handle=fig, data_label='$n=16$', line_style='bo-')
fdsplotlib.plot_to_fig(x_data=r_32, y_data=Y_superbee_32, figure_handle=fig, data_label='$n=32$', line_style='m*-')
fdsplotlib.plot_to_fig(x_data=r_64, y_data=Y_superbee_64, figure_handle=fig, data_label='$n=64$', line_style='r^-')
fdsplotlib.plot_to_fig(x_data=r_128, y_data=Y_superbee_128, figure_handle=fig, data_label='$n=128$', line_style='gsq-')

plt.savefig(os.path.join(pltdir, 'soborot_superbee_square_wave.pdf'), format='pdf')
plt.close()

# SQUARE WAVE: GODUNOV

Y_godunov_16 = read_tracers(os.path.join(outdir, 'soborot_godunov_square_wave_16_devc.csv'), 16)
Y_godunov_32 = read_tracers(os.path.join(outdir, 'soborot_godunov_square_wave_32_devc.csv'), 32)
Y_godunov_64 = read_tracers(os.path.join(outdir, 'soborot_godunov_square_wave_64_devc.csv'), 64)
Y_godunov_128 = read_tracers(os.path.join(outdir, 'soborot_godunov_square_wave_128_devc.csv'), 128)

fig = fdsplotlib.plot_to_fig(x_data=r_exact, y_data=Y_exact, data_label='Exact', line_style='k-',
                  x_label='Radial Position (m)', y_label='Scalar Mass Fraction',
                  plot_title='Godunov',
                  x_min=0, x_max=1, y_min=0, y_max=1.2)
fdsplotlib.plot_to_fig(x_data=r_16, y_data=Y_godunov_16, figure_handle=fig, data_label='$n=16$', line_style='bo-')
fdsplotlib.plot_to_fig(x_data=r_32, y_data=Y_godunov_32, figure_handle=fig, data_label='$n=32$', line_style='m*-')
fdsplotlib.plot_to_fig(x_data=r_64, y_data=Y_godunov_64, figure_handle=fig, data_label='$n=64$', line_style='r^-')
fdsplotlib.plot_to_fig(x_data=r_128, y_data=Y_godunov_128, figure_handle=fig, data_label='$n=128$', line_style='gsq-')

plt.savefig(os.path.join(pltdir, 'soborot_godunov_square_wave.pdf'), format='pdf')
plt.close()

# Square wave error calculations

Y_exact_16 = np.interp(r_16, r_exact, Y_exact)
Y_exact_32 = np.interp(r_32, r_exact, Y_exact)
Y_exact_64 = np.interp(r_64, r_exact, Y_exact)
Y_exact_128 = np.interp(r_128, r_exact, Y_exact)

e_charm_16 = np.linalg.norm(Y_charm_16 - Y_exact_16, 1) / len(r_16)
e_charm_32 = np.linalg.norm(Y_charm_32 - Y_exact_32, 1) / len(r_32)
e_charm_64 = np.linalg.norm(Y_charm_64 - Y_exact_64, 1) / len(r_64)
e_charm_128 = np.linalg.norm(Y_charm_128 - Y_exact_128, 1) / len(r_128)

e_superbee_16 = np.linalg.norm(Y_superbee_16 - Y_exact_16, 1) / len(r_16)
e_superbee_32 = np.linalg.norm(Y_superbee_32 - Y_exact_32, 1) / len(r_32)
e_superbee_64 = np.linalg.norm(Y_superbee_64 - Y_exact_64, 1) / len(r_64)
e_superbee_128 = np.linalg.norm(Y_superbee_128 - Y_exact_128, 1) / len(r_128)
e_superbee_128_1mesh = np.linalg.norm(Y_superbee_128_1mesh - Y_exact_128, 1) / len(r_128)

if abs(e_superbee_128_1mesh - e_superbee_128) > 1.e-10:
    print('Error: soborot_superbee_square_wave_128 single mesh and multi-mesh out of tolerance')

e_godunov_16 = np.linalg.norm(Y_godunov_16 - Y_exact_16, 1) / len(r_16)
e_godunov_32 = np.linalg.norm(Y_godunov_32 - Y_exact_32, 1) / len(r_32)
e_godunov_64 = np.linalg.norm(Y_godunov_64 - Y_exact_64, 1) / len(r_64)
e_godunov_128 = np.linalg.norm(Y_godunov_128 - Y_exact_128, 1) / len(r_128)

dx = L / np.array([16.0, 32.0, 64.0, 128.0])

fig = fdsplotlib.plot_to_fig(x_data=dx, y_data=dx**0.5, data_label=r'$O(\delta x^{1/2})$', line_style='k-.',
                  x_label='Grid Spacing (m)', y_label='L$_2$ Error', plot_type='loglog',
                  plot_title='Square Wave Error',
                  x_min=np.min(dx), x_max=0.1, y_min=np.min(dx), y_max=np.max(dx**0.5))
fdsplotlib.plot_to_fig(x_data=dx, y_data=dx, figure_handle=fig, data_label=r'$O(\delta x)$', line_style='k--')
fdsplotlib.plot_to_fig(x_data=dx, y_data=np.array([e_godunov_16, e_godunov_32, e_godunov_64, e_godunov_128]),
                  figure_handle=fig, data_label='Godunov', line_style='ksq-')
fdsplotlib.plot_to_fig(x_data=dx, y_data=np.array([e_superbee_16, e_superbee_32, e_superbee_64, e_superbee_128]),
                  figure_handle=fig, data_label='Superbee', line_style='k*-')
fdsplotlib.plot_to_fig(x_data=dx, y_data=np.array([e_charm_16, e_charm_32, e_charm_64, e_charm_128]),
                  figure_handle=fig, data_label='CHARM', line_style='ko-')

if e_godunov_128 > 4.2e-02:
    print('Error: soborot_godunov_square_wave_128 out of tolerance')
if e_charm_128 > 1.4e-02:
    print('Error: soborot_charm_square_wave_128 out of tolerance')
if e_superbee_128 > 1.0e-02:
    print('Error: soborot_superbee_square_wave_128 out of tolerance')

plt.savefig(os.path.join(pltdir, 'soborot_square_wave_error.pdf'), format='pdf')
plt.close()

# COSINE WAVE: CHARM

Y_charm_16 = read_tracers(os.path.join(outdir, 'soborot_charm_cos_wave_16_devc.csv'), 16)
Y_charm_32 = read_tracers(os.path.join(outdir, 'soborot_charm_cos_wave_32_devc.csv'), 32)
Y_charm_64 = read_tracers(os.path.join(outdir, 'soborot_charm_cos_wave_64_devc.csv'), 64)
Y_charm_128 = read_tracers(os.path.join(outdir, 'soborot_charm_cos_wave_128_devc.csv'), 128)

Y_exact_cos = np.zeros_like(r_exact)
i_range = np.where((r_exact > 0.25) & (r_exact < 0.75))[0]
Y_exact_cos[i_range] = 0.5 * (1.0 + np.cos(4.0 * np.pi * r_exact[i_range]))

fig = fdsplotlib.plot_to_fig(x_data=r_exact, y_data=Y_exact_cos, data_label='Exact', line_style='k-',
                  x_label='Radial Position (m)', y_label='Scalar Mass Fraction', plot_title='CHARM',
                  x_min=0, x_max=1, y_min=0, y_max=1.2)
fdsplotlib.plot_to_fig(x_data=r_16, y_data=Y_charm_16, figure_handle=fig, data_label='$n=16$', line_style='bo-')
fdsplotlib.plot_to_fig(x_data=r_32, y_data=Y_charm_32, figure_handle=fig, data_label='$n=32$', line_style='m*-')
fdsplotlib.plot_to_fig(x_data=r_64, y_data=Y_charm_64, figure_handle=fig, data_label='$n=64$', line_style='r^-')
fdsplotlib.plot_to_fig(x_data=r_128, y_data=Y_charm_128, figure_handle=fig, data_label='$n=128$', line_style='gsq-')

plt.savefig(os.path.join(pltdir, 'soborot_charm_cos_wave.pdf'), format='pdf')
plt.close()

# COSINE WAVE: SUPERBEE 

Y_superbee_16 = read_tracers(os.path.join(outdir, 'soborot_superbee_cos_wave_16_devc.csv'), 16)
Y_superbee_32 = read_tracers(os.path.join(outdir, 'soborot_superbee_cos_wave_32_devc.csv'), 32)
Y_superbee_64 = read_tracers(os.path.join(outdir, 'soborot_superbee_cos_wave_64_devc.csv'), 64)
Y_superbee_128 = read_tracers(os.path.join(outdir, 'soborot_superbee_cos_wave_128_devc.csv'), 128)

fig = fdsplotlib.plot_to_fig(x_data=r_exact, y_data=Y_exact_cos, data_label='Exact', line_style='k-',
                  x_label='Radial Position (m)', y_label='Scalar Mass Fraction', plot_title='Superbee',
                  x_min=0, x_max=1, y_min=0, y_max=1.2)
fdsplotlib.plot_to_fig(x_data=r_16, y_data=Y_superbee_16, figure_handle=fig, data_label='$n=16$', line_style='bo-')
fdsplotlib.plot_to_fig(x_data=r_32, y_data=Y_superbee_32, figure_handle=fig, data_label='$n=32$', line_style='m*-')
fdsplotlib.plot_to_fig(x_data=r_64, y_data=Y_superbee_64, figure_handle=fig, data_label='$n=64$', line_style='r^-')
fdsplotlib.plot_to_fig(x_data=r_128, y_data=Y_superbee_128, figure_handle=fig, data_label='$n=128$', line_style='gsq-')

plt.savefig(os.path.join(pltdir, 'soborot_superbee_cos_wave.pdf'), format='pdf')
plt.close()

# Cosine wave error calculations

Y_exact_16 = np.interp(r_16, r_exact, Y_exact_cos)
Y_exact_32 = np.interp(r_32, r_exact, Y_exact_cos)
Y_exact_64 = np.interp(r_64, r_exact, Y_exact_cos)
Y_exact_128 = np.interp(r_128, r_exact, Y_exact_cos)

e_charm_16 = np.linalg.norm(Y_charm_16 - Y_exact_16, 1) / len(r_16)
e_charm_32 = np.linalg.norm(Y_charm_32 - Y_exact_32, 1) / len(r_32)
e_charm_64 = np.linalg.norm(Y_charm_64 - Y_exact_64, 1) / len(r_64)
e_charm_128 = np.linalg.norm(Y_charm_128 - Y_exact_128, 1) / len(r_128)

e_superbee_16 = np.linalg.norm(Y_superbee_16 - Y_exact_16, 1) / len(r_16)
e_superbee_32 = np.linalg.norm(Y_superbee_32 - Y_exact_32, 1) / len(r_32)
e_superbee_64 = np.linalg.norm(Y_superbee_64 - Y_exact_64, 1) / len(r_64)
e_superbee_128 = np.linalg.norm(Y_superbee_128 - Y_exact_128, 1) / len(r_128)

fig = fdsplotlib.plot_to_fig(x_data=dx, y_data=dx, data_label=r'$O(\delta x)$', line_style='k-.',
                  plot_title='Cosine Wave Error',
                  x_label='Grid Spacing (m)', y_label='L$_2$ Error', plot_type='loglog',
                  x_min=np.min(dx), x_max=0.1, y_min=np.min(dx**2), y_max=np.max(dx))
fdsplotlib.plot_to_fig(x_data=dx, y_data=dx**2, figure_handle=fig, data_label=r'$O(\delta x^2)$', line_style='k--')
fdsplotlib.plot_to_fig(x_data=dx, y_data=np.array([e_charm_16, e_charm_32, e_charm_64, e_charm_128]),
                  figure_handle=fig, data_label='CHARM', line_style='ko-')
fdsplotlib.plot_to_fig(x_data=dx, y_data=np.array([e_superbee_16, e_superbee_32, e_superbee_64, e_superbee_128]),
                  figure_handle=fig, data_label='Superbee', line_style='k*-')

if e_charm_128 > 2.9e-04:
    print('Error: soborot_charm_cos_wave_128 out of tolerance')
if e_superbee_128 > 6.4e-04:
    print('Error: soborot_superbee_cos_wave_128 out of tolerance')

plt.savefig(os.path.join(pltdir, 'soborot_cos_wave_error.pdf'), format='pdf')
plt.close()

