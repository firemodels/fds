#!/bin/bash
export QFDS=./get_singlecasetime_fds.sh
../FDS_Cases.sh | sort -k 1 -n
