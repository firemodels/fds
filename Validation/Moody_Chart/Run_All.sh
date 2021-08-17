#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N16_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N32_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N64_mu0125.fds
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N64_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N8_mu025.fds

$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p1_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p1_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p0001_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p001_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p01_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p1_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p1_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p02_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR s=p02_dpdx=-100_N16.fds

$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_a_10.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_b_10.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_c_10.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_d_10.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_e_10.fds

$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_a_20.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_b_20.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_c_20.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_d_20.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_pressure_drop_e_20.fds

$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_fire_10.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_fire_20.fds
$QFDS -p 80 $DEBUG $QUEUE -d $INDIR tunnel_fire_40.fds

$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-0p01_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-0p01_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-0p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-100_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-1_N32.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_moody_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_poiseuille_N16_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_poiseuille_N32_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_poiseuille_N64_mu0125.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_poiseuille_N64_mu025.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_poiseuille_N8_mu025.fds

$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p1_dpdx=-1_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p1_dpdx=-100_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-100_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p0001_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p001_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p01_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p1_dpdx=-p0001_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p1_dpdx=-p01_N8.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p02_dpdx=-1_N16.fds
$QFDS $DEBUG $QUEUE -d $INDIR geom_s=p02_dpdx=-100_N16.fds

echo FDS cases submitted







