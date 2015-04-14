#!/bin/bash
export QFDS=./get_casetime.sh
./FDS_Cases.sh | sort -k 1 -n
