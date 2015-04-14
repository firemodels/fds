#!/bin/bash
export QFDS=./get_singlecasetime.sh
export RUNCFAST=./get_singlecasetime.sh
export RUNTFDS=./get_singlecasetime.sh
./SMV_Cases.sh | sort -k 1 -n
