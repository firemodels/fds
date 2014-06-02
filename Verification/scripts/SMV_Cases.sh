#!/bin/bash

$RUNFDS Immersed_Boundary_Method geom_arch
$RUNFDS Immersed_Boundary_Method geom_azim
$RUNFDS Immersed_Boundary_Method geom_elev
$RUNFDS Immersed_Boundary_Method geom_obst
$RUNFDS Immersed_Boundary_Method geom_scale
$RUNFDS Immersed_Boundary_Method geom_simple
$RUNFDS Immersed_Boundary_Method geom_sphere1a
$RUNFDS Immersed_Boundary_Method geom_sphere1b
$RUNFDS Immersed_Boundary_Method geom_sphere1c
$RUNFDS Immersed_Boundary_Method geom_sphere1d
$RUNFDS Immersed_Boundary_Method geom_sphere1e
$RUNFDS Immersed_Boundary_Method geom_sphere1f
$RUNFDS Immersed_Boundary_Method geom_sphere2
$RUNFDS Immersed_Boundary_Method geom_sphere3a
$RUNFDS Immersed_Boundary_Method geom_sphere3b
$RUNFDS Immersed_Boundary_Method geom_sphere3c
$RUNFDS Immersed_Boundary_Method geom_sphere3d
$RUNFDS Immersed_Boundary_Method geom_sphere3e
$RUNFDS Immersed_Boundary_Method geom_sphere3f
$RUNFDS Immersed_Boundary_Method geom_sphere_fire
$RUNFDS Immersed_Boundary_Method geom_terrain
$RUNFDS Immersed_Boundary_Method geom_texture
$RUNFDS Immersed_Boundary_Method geom_texture2
$RUNFDS Immersed_Boundary_Method geom_texture3a
$RUNFDS Immersed_Boundary_Method geom_texture3b
$RUNFDS Immersed_Boundary_Method geom_time
$RUNFDS Immersed_Boundary_Method geom_time2
$RUNFDS Immersed_Boundary_Method geom_time3
$RUNFDS Immersed_Boundary_Method geom_time4
$RUNFDS Immersed_Boundary_Method geom_time5

$RUNFDS Visualization cell_test
$RUNCFAST Visualization cfast_test
$RUNFDS Visualization colorbar
$RUNFDS Visualization colorconv
$RUNFDS Visualization fed_test
$RUNFDS Visualization mplume5c8
$RUNFDS Visualization objects_dynamic
$RUNFDS Visualization objects_elem
$RUNFDS Visualization objects_static
$RUNFDS Visualization plume5c
$RUNFDS Visualization plume5cdelta
$RUNFDS Visualization plumeiso
$RUNFDS Visualization plume5c_bounddef
$RUNFDS Visualization script_test
$RUNFDS Visualization script_slice_test
$RUNFDS Visualization sillytexture
$RUNFDS Visualization slicemask
$RUNFDS Visualization smoke_sensor
$RUNFDS Visualization smoke_test
$RUNFDS Visualization smoke_test2
$RUNFDS Visualization sprinkler_many
$RUNFDS Visualization thouse5
$RUNFDS Visualization thouse5delta
$RUNFDS Visualization transparency
$RUNFDS Visualization vcirctest
$RUNFDS Visualization vcirctest2
$RUNTFDS Visualization version

$RUNFDS WUI levelset1
$RUNWFDS WUI tree_one
#$RUNFDS WUI fire_line
$RUNFDS WUI BT10m_2x2km_LS 1
$RUNFDS WUI wind_test1
$RUNFDS WUI tree_test2
$RUNFDS WUI hill_structure
#$RUNFDS WUI pine_tree

# $RUNFDS WUI onetree_surf_1mesh
