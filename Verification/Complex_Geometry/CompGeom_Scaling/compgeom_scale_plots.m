% compgeom_scale_plots:
%
% -------------------------------------------------------------------------
close all
clear all
clc



fdsdir  = '/Volumes/mnv/FIREMODELS_FORK/fds/';
casedir = 'Verification/Complex_Geometry/CompGeom_Scaling/';
basedir = [fdsdir casedir];
addpath(basedir)

% Constants:
N_CUTCELLS                =  1;
N_INB_CUTFACES            =  2;
N_REG_CUTFACES            =  3; 
SET_CUTCELLS              =  4;
GET_BODINT_PLANE          =  5;
GET_X2_INTERSECTIONS      =  6;
GET_X2_VERTVAR            =  7;
GET_CARTEDGE_CUTEDGES     =  8;
GET_BODX2X3_INTERSECTIONS =  9;
GET_CARTFACE_CUTEDGES     = 10; 
GET_CARTCELL_CUTEDGES     = 11; 
GET_CARTFACE_CUTFACES     = 12;
GET_CARTCELL_CUTFACES     = 13;
GET_CARTCELL_CUTCELLS     = 14;


% Scaling files for increasing number of triangles, serial run:
N_TRIANG = [ 320 1280 5120 20480 81920 327680 ]'; 
serial_files_t = {'compgeom_scale_64x64_1mesh_320T_cc_cpu_000.csv'  ; ...
                  'compgeom_scale_64x64_1mesh_1280T_cc_cpu_000.csv' ; ...
                  'compgeom_scale_64x64_1mesh_5120T_cc_cpu_000.csv' ; ...
                  'compgeom_scale_64x64_1mesh_20480T_cc_cpu_000.csv'; ...
                  'compgeom_scale_64x64_1mesh_81920T_cc_cpu_000.csv'; ...
                  'compgeom_scale_64x64_1mesh_327680T_cc_cpu_000.csv'};
                  
serial_tags_t  = {'64^3_1mesh_320T'  ; ...
                  '64^3_1mesh_1280T' ; ...
                  '64^3_1mesh_5120T' ; ...
                  '64^3_1mesh_20480T'; ...
                  '64^3_1mesh_81920T'; ...
                  '64^3_1mesh_327680T'};
              
n_files_t = length(serial_files_t);

for ifile=1:n_files_t
    filename=serial_files_t{ifile};
    TIME_T(ifile).filename = filename;
    TIME_T(ifile).tagname  = serial_tags_t{ifile};
    
    % Load number of geometric entities and timings:
    [fid]=fopen([basedir filename],'r');
    
    line =fgetl(fid);
    line2=fgetl(fid)
    vec =str2num(line2)

    TIME_T(ifile).geomelem = vec(N_CUTCELLS:N_REG_CUTFACES);
    TIME_T(ifile).time     = vec(SET_CUTCELLS:GET_CARTCELL_CUTCELLS);
    TIME_T2(ifile,N_CUTCELLS:N_REG_CUTFACES) = TIME_T(ifile).geomelem;
    TIME_T2(ifile,SET_CUTCELLS:GET_CARTCELL_CUTCELLS) = TIME_T(ifile).time;
    fclose(fid);
end

figure
hold on
plot(N_TRIANG(:),TIME_T2(:,GET_BODINT_PLANE))
%plot(N_TRIANG,TIME_T2(:,GET_CARTCELL_CUTCELLS))
plot(N_TRIANG,0.00001*(N_TRIANG))
plot(N_TRIANG,0.1*sqrt(N_TRIANG))
xlabel('N triangles')
ylabel('Wall time [sec]')
set(gca,'Xscale','log','Yscale','log')
box on


