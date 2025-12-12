
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/NIST_Pool_Fires/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/NIST_Pool_Fires/'

chids = ['NIST_Acetone_Prescribed_0p5cm',
         'NIST_Ethanol_Prescribed_0p5cm',
         'NIST_Heptane_Prescribed_0p5cm',
         'NIST_Methane_Prescribed_0p5cm',
         'NIST_Methanol_Prescribed_0p5cm',
         'NIST_Propane_20kW_Prescribed_0p5cm',
         'NIST_Propane_34kW_Prescribed_0p5cm',
         'NIST_Propane_50kW_Prescribed_0p5cm']

fuel = ['Acetone','Ethanol','Heptane','Methane','Methanol','Propane 20 kW','Propane 34 kW','Propane 50 kW']
label = ['Acetone','Ethanol','Heptane','Methane','Methanol','Propane_20','Propane_34','Propane_50']
ideal = [2.45,2.41,2.8,2.48,2.49,2.22,2.39,2.39] # Heptane value based on theory, not measurement.
i = -1

for chid in chids:

    i = i + 1
    git_file = outdir + chid + '_git.txt'
    version_string = fdsplotlib.get_version_string(git_file)
    
    # Read CSV file, skipping the first two rows
    df = pd.read_csv(outdir + chid + '_devc.csv', skiprows=20000, header=None)
    
    # Extract time and data columns
    time = df.iloc[:, 0].values
    data = df.iloc[:, 1].values
    
    # Calculate the sampling rate
    dt = np.mean(np.diff(time))  # Average time step
    sampling_rate = 1.0 / dt
    
    # Compute FFT
    n = len(data)
    fft_values = np.fft.fft(data)
    fft_freq = np.fft.fftfreq(n, dt)
    
    # Get positive frequencies only
    positive_freq_idx = fft_freq > 0
    frequencies = fft_freq[positive_freq_idx]
    magnitude = np.abs(fft_values[positive_freq_idx])
    
    # Create the frequency spectrum plot using plot_to_fig
    fig = fdsplotlib.plot_to_fig(frequencies, magnitude, marker_style='k-',
                                 x_label='Frequency (Hz)', y_label='Magnitude',
                                 x_min=0, x_max=6, y_min=0,
                                 revision_label=version_string,
                                 plot_title=fuel[i])
    fdsplotlib.plot_to_fig([ideal[i],ideal[i]],[0,100000], marker_style='k--', figure_handle=fig)
    
    if i==2:
       ax = plt.gca()
       ax.text(0.3, 3000, 'Not measured; theory only')

    fig.savefig(pltdir + label[i] + '_frequency_spectrum.pdf', format='pdf')

