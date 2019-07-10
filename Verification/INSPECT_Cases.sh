#!/bin/bash

# cases in this script are run with the openmp thread checker version of fds
# these case are runot run in debug mode (as are other VV cases)

$QFDS -o 4 -d Thread_Check -x results race_test_4.fds
