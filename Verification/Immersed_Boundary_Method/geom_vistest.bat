@echo off
:: set smv=call ..\..\..\bin\smv_r -runscript -bindir ..\..\SMV\for_bundle
set smv=..\..\SMV\Build\intel_win_64\smokeview_win_64 -runscript -bindir ..\..\SMV\for_bundle
%smv% geom_arch
%smv% geom_azim
%smv% geom_elev
%smv% geom_obst
%smv% geom_scale
%smv% geom_simple
%smv% geom_sphere
%smv% geom_sphere2
%smv% geom_terrain
%smv% geom_texture
%smv% geom_texture2
%smv% geom_texture3a
%smv% geom_texture3b
%smv% geom_time
%smv% geom_time2
%smv% geom_time3
%smv% geom_time4
%smv% geom_time5
