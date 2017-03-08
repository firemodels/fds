# K. Overholt
# 9-3-2013
# validation_uncertainty.py

# This script generates input files for the quantification of uncertainty of
# the THIEF model (CAROLFIRE cases) using the correlations program

from __future__ import division
import numpy as np
import pandas as pd

#  =================
#  = User settings =
#  =================

# mode = 'write' input file or 'read' results
mode = 'write'

# Number of MC iterations
n_samples = 10000

# Cable failure temperature (C)
failure_temp = 400

# Input and output file locations
input_file = 'validation_uncertainty_THIEF_inputs.txt'
output_file = 'validation_uncertainty_THIEF_outputs'

#  ====================
#  = Write input file =
#  ====================

if mode == 'write':
    # Input parameter specifications/distributions
    D = np.random.uniform(16.25, 16.35, n_samples)
    MASS_PER_LENGTH = 0.529
    JACKET_THICKNESS = np.random.uniform(1.45, 1.55, n_samples)
    TMP_A = 20
    TMP_RAMP = TMP_A + np.random.normal(450-TMP_A, 0.025*(450-TMP_A), n_samples)
    T_END = 1800

    # Generate correlation input file
    f = open(input_file, 'w')
    for i in range(len(TMP_RAMP)):
        output_string = "{} D={:.2f}, MASS_PER_LENGTH={:.3f}, JACKET_THICKNESS={:.2f}, TMP_A={:.0f}, T_RAMP(0:5)=0,70,120,180,240,300, TMP_RAMP(0:5)={:.1f},{:.1f},{:.1f},{:.1f},{:.1f},{:.1f}, T_END={:.0f}, OUTPUT_FILE='{}' /\n".format('&THIEF', D[i], MASS_PER_LENGTH, JACKET_THICKNESS[i], TMP_A, TMP_A, TMP_RAMP[i], TMP_RAMP[i], TMP_RAMP[i], TMP_RAMP[i],TMP_RAMP[i], T_END, output_file + '_' + str(i) + '.csv')
        f.write(output_string)
    f.close()

#  ================
#  = Read results =
#  ================

if mode == 'read':
    activation_times = np.array([])

    for i in range(n_samples):
        output_file = 'validation_uncertainty_THIEF_outputs_' + str(i) + '.csv'
        data = pd.read_table(output_file, sep=',')

        for j in range(len(data)):
            if (data['Cable Temp'][j] >= failure_temp):
                activation_times = np.append(activation_times, data['Time'][j])
                break

    print 'Times:', activation_times
    print 'Standard Deviation:', activation_times.std()
    print 'Mean:', activation_times.mean()
    print 'Sigma_E:', (activation_times.std() / activation_times.mean())
