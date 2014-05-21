#!/bin/bash
#  1.  after making changes to the template input file in
#      makecase.sh, run this script
#  2.  then commit the updated test...fds input files to the repository

./makecase.sh 64 test64a
./makecase.sh 64 test64b
./makecase.sh 64 test64c
./makecase.sh 64 test64d
./makecase.sh 64 test64e
./makecase.sh 64 test64f
./makecase.sh 64 test64g
./makecase.sh 64 test64h
./makecase.sh 128 test128a
./makecase.sh 128 test128b
./makecase.sh 128 test128c
./makecase.sh 128 test128d
./makecase.sh 128 test128e
./makecase.sh 128 test128f
./makecase.sh 128 test128g
./makecase.sh 128 test128h
