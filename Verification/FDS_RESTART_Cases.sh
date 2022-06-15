#!/bin/bash

$QFDS -p 4 -d Restart device_restart_b.fds
$QFDS -d Restart restart_test1b.fds
$QFDS -d Restart pyro3d_restart_b.fds
$QFDS -d Restart geom_restart_b.fds
$QFDS -d Restart clocks_restart_b.fds
