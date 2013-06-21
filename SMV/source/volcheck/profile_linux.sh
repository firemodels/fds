#!/bin/bash
./make_linuxp.sh
./volcheck_linux_64p
gprof ./volcheck_linux_64p > results.out
