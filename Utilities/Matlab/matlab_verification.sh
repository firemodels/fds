#!/bin/bash

# Replace LaTeX with TeX for Interpreter in plot_style.m
# This allows displayless automatic Matlab plotting
# Otherwise Matlab crashes due to a known bug
cd scripts
sed -i 's/LaTeX/TeX/g' plot_style.m 
cd ..

matlab -r "try, disp('Running Matlab Verification script'), FDS_verification_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"

# Restore LaTeX as plot_style interpreter
cd scripts
sed -i 's/TeX/LaTeX/g' plot_style.m
cd ..
