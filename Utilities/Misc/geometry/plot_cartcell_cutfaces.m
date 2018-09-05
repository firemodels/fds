% plot_cartcell_cutfaces:
% Routine to visualize problematic INBOUNDARY cut-face definition in a
% particular Cartesian cell:
%
% -------------------------------------------------------------------------
close all
clear all
clc


IAXIS = 1; JAXIS = 2; KAXIS = 3;
NOD1  = 1; NOD2  = 2; NOD3  = 3;

basedir='/Users/mnv/Documents/FIREMODELS_FORK/fds/test/';
file='Cartcell_cutfaces.dat';

% Load cartesian cell cut-faces intersection data:
[fid]=fopen([basedir file],'r');
% I,J,K: 
[cline]=fgetl(fid);
vec    = str2num(fgetl(fid));
I = vec(1); 
J = vec(2);
K = vec(3);
% XC,YC,ZC and DXC,DYC,DZC:
[cline]=fgetl(fid);
for x1axis=IAXIS:KAXIS
    vec    = str2num(fgetl(fid));
    xyzc(x1axis)  = vec(1);
    dxyzc(x1axis) = vec(2);
end
% NVERT,NSEG,NSEG_FACE,COUNTR,NSEG_LEFT:
[cline]=fgetl(fid);
vec       = str2num(fgetl(fid));
NVERT     = vec(1);
NSEG      = vec(2);
NSEG_FACE = vec(3);
COUNTR    = vec(4);
NSEG_LEFT = vec(5);
% XYZVERT:
[cline]=fgetl(fid);
for ivert = 1:NVERT
    vec   = str2num(fgetl(fid));
    XYZVERT(ivert,IAXIS:KAXIS) = vec(2:4);
end
% SEG_CELL:
[cline]=fgetl(fid);
for iseg  = 1:NSEG
    vec   = str2num(fgetl(fid));
    SEG_CELL(iseg,1:6) = vec(2:7);
end
% Original SEG_FACE:
[cline]=fgetl(fid);
for iseg  = 1:NSEG_FACE
    vec   = str2num(fgetl(fid));
    SEG_FACE(iseg,NOD1:NOD2) = vec(2:3);
end
% SEG_FACE2 computed so far:
[cline]=fgetl(fid);
for iseg  = 1:COUNTR
    vec   = str2num(fgetl(fid));
    SEG_FACE2(iseg,NOD1:NOD2) = vec(2:3);
end
% 'ICF,BOD_TRI:'
[cline] = fgetl(fid);
vec     = str2num(fgetl(fid));
ICF     = vec(1);
NBODTRI = vec(2);
for idum=1:NBODTRI
    vec     = str2num(fgetl(fid));
    BOD_TRI(idum,1:2) = vec(1:2);
end
fclose(fid);

% Now make the plot:
% Load Geometry:
geom_vert_file='GEOMETRY_0001_VERTS.dat';
geom_face_file='GEOMETRY_0001_FACES.dat';
XYZ=load([basedir geom_vert_file]);
WSELEM=load([basedir geom_face_file]);

figure
hold on
% Geometry:
[hg]=trisurf(WSELEM,XYZ(:,IAXIS),XYZ(:,JAXIS),XYZ(:,KAXIS));
xlabel('X')
ylabel('Y')
zlabel('Z')
box on
grid on
axis equal
view([45 45])


% Plot Cell:
ct=0;
for k=1:2
    for j=1:2
        for i=1:2
            ct=ct+1;
            xyz_cell(ct,IAXIS:KAXIS)=[xyzc(IAXIS)+(i-3/2)*dxyzc(IAXIS) ...
                                      xyzc(JAXIS)+(j-3/2)*dxyzc(JAXIS) ...
                                      xyzc(KAXIS)+(k-3/2)*dxyzc(KAXIS)];
        end
    end
end

seg_list = [ 1 2; 3 4; 5 6; 7 8; 1 3; 2 4; 5 7; 6 8; 1 5; 2 6; 3 7; 4 8];
for iseg  = 1:12
    xyz1 = xyz_cell(seg_list(iseg,NOD1),IAXIS:KAXIS);
    xyz2 = xyz_cell(seg_list(iseg,NOD2),IAXIS:KAXIS);
    plot3([xyz1(IAXIS) xyz2(IAXIS)], ...
          [xyz1(JAXIS) xyz2(JAXIS)], ...
          [xyz1(KAXIS) xyz2(KAXIS)], 'k', 'Linewidth', 1);
end

% Now plot on top the Segments, from SEG_CELL and SEG_FACE
for iseg  = 1:NSEG
    xyz1 = XYZVERT(SEG_CELL(iseg,NOD1),IAXIS:KAXIS);
    xyz2 = XYZVERT(SEG_CELL(iseg,NOD2),IAXIS:KAXIS);
    plot3([xyz1(IAXIS) xyz2(IAXIS)], ...
          [xyz1(JAXIS) xyz2(JAXIS)], ...
          [xyz1(KAXIS) xyz2(KAXIS)], 'b', 'Linewidth', 4);
end
% Original SEG_FACE:
for iseg  = 1:NSEG_FACE
    xyz1 = XYZVERT(SEG_FACE(iseg,NOD1),IAXIS:KAXIS);
    xyz2 = XYZVERT(SEG_FACE(iseg,NOD2),IAXIS:KAXIS);
    plot3([xyz1(IAXIS) xyz2(IAXIS)], ...
          [xyz1(JAXIS) xyz2(JAXIS)], ...
          [xyz1(KAXIS) xyz2(KAXIS)], 'r', 'Linewidth', 4);
end

pause

axis([min(xyz_cell(:,IAXIS))-dxyzc(IAXIS) max(xyz_cell(:,IAXIS))+dxyzc(IAXIS) ...
      min(xyz_cell(:,JAXIS))-dxyzc(JAXIS) max(xyz_cell(:,JAXIS))+dxyzc(JAXIS) ...
      min(xyz_cell(:,KAXIS))-dxyzc(KAXIS) max(xyz_cell(:,KAXIS))+dxyzc(KAXIS)])









