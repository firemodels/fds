#!/bin/bash
#  1.  after making changes to the template input file in
#      makecase.sh, run this script
#  2.  then commit the updated openmp_testxxx.fds input files to the repository

./makecase.sh 64 openmp_test64a
./makecase.sh 64 openmp_test64b
./makecase.sh 64 openmp_test64c
./makecase.sh 64 openmp_test64d
./makecase.sh 64 openmp_test64e
./makecase.sh 64 openmp_test64f
./makecase.sh 64 openmp_test64g
./makecase.sh 64 openmp_test64h
./makecase.sh 128 openmp_test128a
./makecase.sh 128 openmp_test128b
./makecase.sh 128 openmp_test128c
./makecase.sh 128 openmp_test128d
./makecase.sh 128 openmp_test128e
./makecase.sh 128 openmp_test128f
./makecase.sh 128 openmp_test128g
./makecase.sh 128 openmp_test128h
