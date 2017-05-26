#!/bin/bash
#  1.  after making changes to the template input file in
#      makecase.sh, run this script
#  2.  then commit the updated openmp_testxxx.fds input files to the repository

./makecase64.sh 64 openmp_test64a
./makecase64.sh 64 openmp_test64b
./makecase64.sh 64 openmp_test64c
./makecase64.sh 64 openmp_test64d
./makecase64.sh 64 openmp_test64e
./makecase64.sh 64 openmp_test64f
./makecase64.sh 64 openmp_test64g
./makecase64.sh 64 openmp_test64h
./makecase128.sh 128 openmp_test128a
./makecase128.sh 128 openmp_test128b
./makecase128.sh 128 openmp_test128c
./makecase128.sh 128 openmp_test128d
./makecase128.sh 128 openmp_test128e
./makecase128.sh 128 openmp_test128f
./makecase128.sh 128 openmp_test128g
./makecase128.sh 128 openmp_test128h
