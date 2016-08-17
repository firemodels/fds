#!/bin/bash
fds=fds
#fds=../../../../FDS_Compilation/intel_linux_64_db/fds_linux_64_db
#fds=../../../../FDS_Compilation/intel_linux_64/fds_linux_64

$fds blodget.fds

$fds demtest1.fds
$fds demtest2.fds

$fds nist.fds

$fds sugarloaf.fds

$fds tower.fds

$fds trails.fds
$fds trails2.fds
