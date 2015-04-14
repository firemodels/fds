#!/bin/bash
export QFDS=./get_singlecasetime.sh
./FDS_Cases.sh | sort -k 1 -n
