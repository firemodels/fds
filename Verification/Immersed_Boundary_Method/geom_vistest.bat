@echo off

:: repository smokeview
:: set smv=..\..\SMV\Build\intel_win_64\smokeview_win_64 -bindir ..\..\SMV\for_bundle -runscript

:: installed smokeview
set smv=smokeview  -runscript

%smv% geom_arch
%smv% geom_azim
%smv% geom_elev
%smv% geom_obst
%smv% geom_scale
%smv% geom_simple
%smv% geom_sphere1a
%smv% geom_sphere1b
%smv% geom_sphere1c
%smv% geom_sphere1d
%smv% geom_sphere1e
%smv% geom_sphere1f
%smv% geom_sphere2
%smv% geom_sphere3a
%smv% geom_sphere3b
%smv% geom_sphere3c
%smv% geom_sphere3d
%smv% geom_sphere3e
%smv% geom_sphere3f
%smv% geom_sphere_fire
%smv% geom_terrain
%smv% geom_texture
%smv% geom_texture2
%smv% geom_texture3a
%smv% geom_texture3b
%smv% geom_texture4a
%smv% geom_texture4b
%smv% geom_time
%smv% geom_time2
%smv% geom_time3
%smv% geom_time4
%smv% geom_time5
