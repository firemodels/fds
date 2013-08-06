@echo off
matlab -r "try, disp('Running Matlab Validation script'), FDTs_verification_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"