#!/bin/bash

matlab -r "try, disp('Running Matlab Validation script'), FDTs_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"
