% Plot REG-CUTCELL segments for wall model definition:
% Activate #define DEBUG_SET_CUTCELLS
%          #define DEBUG_IBM_INTERPOLATION
% to write out required files.
% -------------------------------------------------------------------------
close all
clear all
clc

IAXIS = 1; JAXIS = 2; KAXIS = 3;
NOD1  = 1; NOD2  = 2; NOD3  = 3;

basedir ='/Users/mnv/Documents/FIREMODELS_FORK/fds/Verification/Complex_Geometry/';
casename='sphere_propane_demo';

N_GEOMETRY=1;

% First load meshes:
[MESHES,NMESHES]=load_meshes(basedir,casename);

newfig=1;
[ierr]=plot_meshes(NMESHES,MESHES,newfig);


% Then load Geometries:
[GEOM,N_GEOMETRY]=load_geometries(basedir,casename);

newfig=0;
[ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig);


% Finally load REG-CUTCELL segments:
[MESHES]=load_rcedges(basedir,casename,NMESHES,MESHES);

[ierr]=plot_rcedges(NMESHES,MESHES,newfig);

return
