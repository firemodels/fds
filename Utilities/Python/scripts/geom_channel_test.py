"""
geom_channel_test.py
Temporary home for checking AREA integration for geom_channel.fds
Original MATLAB script by McDermott (8-11-2021)
"""
import numpy as np
import pandas as pd
import os

def check_geom_channel():
    """Check AREA integration for geom_channel test case"""
    
    fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
    ddir = os.path.join(fds_dir, 'Verification','Complex_Geometry','')    
    devc_file = os.path.join(ddir, 'geom_channel_devc.csv')
    try:
        M = pd.read_csv(devc_file, skiprows=1)  # Skip units row, use column names as header
    except FileNotFoundError:
        print(f'Error: File {devc_file} does not exist. Skipping case.')
        return False
    
    U0 = M['U0'].iloc[-1]
    A0 = M['A0'].iloc[-1]
    U1 = M['U1'].iloc[-1]
    A1 = M['A1'].iloc[-1]
    U2 = M['U2'].iloc[-1]
    A2 = M['A2'].iloc[-1]
    U3 = M['U3'].iloc[-1]
    A3 = M['A3'].iloc[-1]
    
    error_tolerance = 1e-6
    all_passed = True    
    A0_error = abs(A0 - 1)
    if A0_error > error_tolerance:
        print(f'Python Warning: geom_channel.fds A0_error = {A0_error:.6e}')
        all_passed = False
    A1_error = abs(A1 - 1)
    if A1_error > error_tolerance:
        print(f'Python Warning: geom_channel.fds A1_error = {A1_error:.6e}')
        all_passed = False
    A2_error = abs(A2 - 1)
    if A2_error > error_tolerance:
        print(f'Python Warning: geom_channel.fds A2_error = {A2_error:.6e}')
        all_passed = False
    A3_error = abs(A3 - 1)
    if A3_error > error_tolerance:
        print(f'Python Warning: geom_channel.fds A3_error = {A3_error:.6e}')
        all_passed = False
    return all_passed

if __name__ == '__main__':
    """Main execution block"""
    success = check_geom_channel()

