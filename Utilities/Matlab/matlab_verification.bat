@echo off
 matlab -r "try, disp('Running Matlab Verification script'), FDS_verification_script, catch, disp('Matlab error'), err = lasterror, err.message, err.stack, end, exit"