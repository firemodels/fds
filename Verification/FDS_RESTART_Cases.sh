#!/bin/bash

$QFDS -p 4 -d Restart device_restart_b.fds
$QFDS -d Restart restart_test1b.fds
$QFDS -d Restart geom_restart_b.fds
$QFDS -d Restart geom_ls_restart_b.fds
$QFDS -d Restart clocks_restart_b.fds
$QFDS -p 4 -d Restart restart_ulmat_b.fds
$QFDS -p 8 -d Restart csvf_restart_b.fds
