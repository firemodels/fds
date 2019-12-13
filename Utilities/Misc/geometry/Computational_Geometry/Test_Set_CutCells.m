close all 
clear all
clc

% Declarations, call set constants:
global MAX_DIM IAXIS JAXIS KAXIS NOD1 NOD2 NOD3 NOD4
global NMESHES MESHES BODINT_PLANE
global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_BDNUM_CRS_AUX IBM_IS_CRS2_AUX
global NM N_GEOMETRY GEOM IBM_INBOUNDARY

global XCELL DXCELL XFACE DXFACE 
global YCELL DYCELL YFACE DYFACE
global ZCELL DZCELL ZFACE DZFACE


plot_cutedges=false;

[ierr]=SET_CONSTANTS();

% Case name and directory:
addpath ../
basedir = '/Users/mnv/Documents/FIREMODELS_FORK/fds/GEOM_Intersection/';
casename= 'two_spheres';

% Load and plot Geometries:
[GEOM,N_GEOMETRY]=load_geometries(basedir,casename);
newfig =1;
[ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig);

% Load and plot Meshes:
[MESHES,NMESHES]=load_meshes(basedir,casename);
newfig=0;
[ierr]=plot_meshes(NMESHES,MESHES,newfig);

tic
[ierr]=SET_CUTCELLS_3D(basedir,casename,plot_cutedges);
toc

figure
axis equal; box on;
xlabel('X')
ylabel('Y')
zlabel('Z')
NM=1;
for ICF=1:MESHES(NM).N_CUTFACE_MESH
   if(MESHES(NM).CUT_FACE(ICF).STATUS ~= IBM_INBOUNDARY); continue; end
   NFACE=MESHES(NM).CUT_FACE(ICF).NFACE;
   XYZVERT = MESHES(NM).CUT_FACE(ICF).XYZVERT;
   for JCF=1:NFACE
       NELEM  = MESHES(NM).CUT_FACE(ICF).CFELEM(1,JCF);
       CFELEM = MESHES(NM).CUT_FACE(ICF).CFELEM(2:NELEM+1,JCF);
       
       [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'b');
       set(hp,'FaceAlpha',0.5)
       
   end
end

return