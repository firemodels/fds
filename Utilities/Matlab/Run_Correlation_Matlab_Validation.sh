#!/bin/bash

matlab -r "try, disp('Running Matlab Validation script'), Correlation_validation_script, catch, disp('Error'), err = lasterror, err.message, err.stack, end, exit"
