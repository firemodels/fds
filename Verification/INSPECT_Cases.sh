#!/bin/bash

# cases in this script are run with the openmp thread checker version of fds
# these case are runot run in debug mode (as are other VV cases)

$QFDS -d Thread_Check -x Thread_Check/results thread_check.fds
