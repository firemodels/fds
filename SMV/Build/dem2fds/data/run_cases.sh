#!/bin/bash
fds=fds
#fds=../../../../FDS_Compilation/intel_linux_64_db/fds_linux_64_db
#fds=../../../../FDS_Compilation/intel_linux_64/fds_linux_64

case=blodget
fds ${case}.fds

case=nist
fds ${case}.fds

case=sugarloaf
fds ${case}.fds

case=trails
fds ${case}.fds

case=test1x1
$fds ${case}.fds

case=test1x2
$fds ${case}.fds

case=test2x2
$fds ${case}.fds

case=test3x3
$fds ${case}.fds

