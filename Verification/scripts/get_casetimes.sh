#!/bin/bash
export QFDS=./get_casetime.sh
export RUNCFAST=./get_casetime.sh
export RUNTFDS=./get_casetime.sh
./SMV_Cases.sh | sort -k 1 -n
