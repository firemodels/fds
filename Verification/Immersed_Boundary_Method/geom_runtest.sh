#!/bin/bash

# comment or uncomment following lines according to which FDS is running cases

# for local copy in ~/bin directory
# fds=../../../bin/fds6_db

# for repository
# fds=../../FDS_Compilation/intel_osx_64/fds_intel_osx_64
# fds=../../FDS_Compilation/intel_linux_64/fds_intel_linux_64

# for blaze
fds=qfds -r
$fds geom_arch.fds
$fds geom_azim.fds
$fds geom_elev.fds
$fds geom_obst.fds
$fds geom_scale.fds
$fds geom_simple.fds
$fds geom_sphere1a.fds
$fds geom_sphere1b.fds
$fds geom_sphere1c.fds
$fds geom_sphere1d.fds
$fds geom_sphere1e.fds
$fds geom_sphere1f.fds
$fds geom_sphere2.fds
$fds geom_sphere3a.fds
$fds geom_sphere3b.fds
$fds geom_sphere3c.fds
$fds geom_sphere3d.fds
$fds geom_sphere3e.fds
$fds geom_sphere3f.fds
$fds geom_terrain.fds
$fds geom_texture.fds
$fds geom_texture2.fds
$fds geom_texture3a.fds
$fds geom_texture3b.fds
$fds geom_time.fds
$fds geom_time2.fds
$fds geom_time3.fds
$fds geom_time4.fds
$fds geom_time5.fds
$fds geom_volume.fds
