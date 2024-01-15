#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

# FORCED CONVECTION BOUNDARY LAYER SOLUTION

$QFDS $DEBUG -I -p 16 $QUEUE -d $INDIR Pohlhausen_Pr_1.fds
$QFDS $DEBUG -I -p 16 $QUEUE -d $INDIR Pohlhausen_Pr_2.fds
$QFDS $DEBUG -I -p 16 $QUEUE -d $INDIR Pohlhausen_Pr_p5.fds

# FORCED CONVECTION OVER A FLAT PLATE

$QFDS $DEBUG -I $QUEUE -d $INDIR forced_conv_flat_plate_u0p1_25cm.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR forced_conv_flat_plate_u1_25cm.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR forced_conv_flat_plate_u10_25cm.fds

$QFDS $DEBUG -I -p 8 $QUEUE -d $INDIR forced_conv_flat_plate_u0p1_10cm.fds
$QFDS $DEBUG -I -p 8 $QUEUE -d $INDIR forced_conv_flat_plate_u1_10cm.fds
$QFDS $DEBUG -I -p 8 $QUEUE -d $INDIR forced_conv_flat_plate_u10_10cm.fds

$QFDS $DEBUG -I -p 64 $QUEUE -d $INDIR forced_conv_flat_plate_u0p1_2p5cm.fds
$QFDS $DEBUG -I -p 64 $QUEUE -d $INDIR forced_conv_flat_plate_u1_2p5cm.fds
$QFDS $DEBUG -I -p 64 $QUEUE -d $INDIR forced_conv_flat_plate_u10_2p5cm.fds

# NATURAL CONVECTION IN A CLOSED CAVITY

$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_1_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_2_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_3_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_4_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_5_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_6_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_7_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_8_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_9_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_10_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_11_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_12_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_13_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_14_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_15_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_16_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_17_8.fds

$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_1_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_2_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_3_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_4_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_5_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_6_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_7_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_8_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_9_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_10_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_11_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_12_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_13_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_14_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_15_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_16_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_17_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_18_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_19_8.fds

$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_1_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_2_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_3_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_4_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_5_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_6_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_7_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_8_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_9_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_10_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_11_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_12_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_13_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_14_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_15_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_16_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconv_17_16.fds

$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_1_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_2_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_3_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_4_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_5_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_6_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_7_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_8_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_9_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_10_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_11_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_12_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_13_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_14_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_15_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_16_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_17_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_18_16.fds
$QFDS $DEBUG -I -p 16  $QUEUE -d $INDIR natconh_19_16.fds

$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_1_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_2_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_3_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_4_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_5_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_6_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_7_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_8_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_9_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_10_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_11_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_12_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_13_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_14_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_15_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_16_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconv_17_32.fds

$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_1_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_2_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_3_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_4_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_5_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_6_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_7_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_8_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_9_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_10_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_11_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_12_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_13_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_14_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_15_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_16_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_17_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_18_32.fds
$QFDS $DEBUG -I -p 64  $QUEUE -d $INDIR natconh_19_32.fds

# FREE CONVECTION FROM A SPHERE

$QFDS $DEBUG -I -p 8 $QUEUE -d $INDIR free_conv_sphere_1_8.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR free_conv_sphere_2_8.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR free_conv_sphere_3_8.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR free_conv_sphere_4_8.fds

$QFDS $DEBUG -I -p 120 $QUEUE -d $INDIR free_conv_sphere_1_16.fds
$QFDS $DEBUG -I -p 72  $QUEUE -d $INDIR free_conv_sphere_2_16.fds
$QFDS $DEBUG -I -p 72  $QUEUE -d $INDIR free_conv_sphere_3_16.fds
$QFDS $DEBUG -I -p 72  $QUEUE -d $INDIR free_conv_sphere_4_16.fds

# ROTATED NATURAL CONVECTION CASES USING GEOM

$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_1_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_2_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_3_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_4_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_5_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_6_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_7_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_8_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_9_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_10_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_11_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_12_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_13_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_14_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_15_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_16_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_17_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_18_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconh_19_8_rot_18.fds

$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_1_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_2_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_3_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_4_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_5_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_6_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_7_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_8_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_9_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_10_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_11_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_12_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_13_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_14_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_15_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_16_8_rot_18.fds
$QFDS $DEBUG -I -p 2 $QUEUE -d $INDIR natconv_17_8_rot_18.fds

$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_1_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_2_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_3_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_4_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_5_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_6_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_7_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_8_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_9_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_10_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_11_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_12_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_13_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_14_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_15_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_16_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_17_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_18_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconh_19_16_rot_18.fds

$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_1_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_2_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_3_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_4_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_5_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_6_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_7_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_8_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_9_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_10_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_11_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_12_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_13_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_14_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_15_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_16_16_rot_18.fds
$QFDS $DEBUG -I -p 32  $QUEUE -d $INDIR natconv_17_16_rot_18.fds

# NATURAL CONVECTION FROM AN UPWARDS AND DOWNWARDS FACING HOT PLATE

$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_1_16.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_1_32.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_1_8.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_2_16.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_2_32.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_2_8.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_3_16.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_3_32.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_3_8.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_4_16.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_4_32.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_4_8.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_5_16.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_5_32.fds
$QFDS $DEBUG -I -p 8  $QUEUE -d $INDIR nat_conv_hot_plate_5_8.fds

# FORCED CONVECTION FROM AN IMPINGING JET

$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_1e5_coarse.fds
$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_4e5_coarse.fds
$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_1e5_medium.fds
$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_4e5_medium.fds
$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_1e5_fine.fds
$QFDS $DEBUG -I -p 35  $QUEUE -d $INDIR impinging_jet_Re_4e5_fine.fds

echo FDS cases submitted
