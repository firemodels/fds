#!/bin/bash
case=blodget
dem2fds -e $case < ${case}_elevs.csv > ${case}.fds

case=nist
dem2fds -e $case < ${case}_elevs.csv > ${case}.fds

case=sugarloaf
dem2fds -e $case < ${case}_elevs.csv > ${case}.fds

case=test
dem2fds -e $case < ${case}_elevs.csv > ${case}.fds

case=trails
dem2fds -e $case < ${case}_elevs.csv > ${case}.fds
