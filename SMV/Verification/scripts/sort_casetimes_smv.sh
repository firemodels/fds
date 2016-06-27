#!/bin/bash
export QFDS=./get_singlecasetime_smv.sh
export RUNCFAST=./get_singlecasetime_smv.sh
export RUNTFDS=./get_singlecasetime_smv.sh
./SMV_Cases.sh | sort -k 1 -n
