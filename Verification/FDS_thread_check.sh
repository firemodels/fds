#!/bin/bash

# cases in this script are run with the openmp thread checker version of fds
# these case are runot run in debug mode (as are other VV cases)

$QFDS -T -d Thread_Check thread_check.fds
