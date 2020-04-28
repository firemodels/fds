%% ------------------------------------------------------------------------
% Script stl2fds:
%
% This script reads an .stl 3D geometry file, and writes the geometry nodes 
% and surface triangles in FEM_MESH FDS format (ASCII). The result can be
% directly copied to FDS input files.
%
% The geometry must be: A single geometry with water tight surface,
% obtained from a single solid volume on CAD software.
%
% In order to use this script you need to download the free package STL
% File Reader by Eric Johnson, form Matlab Central. The URL is:
%
% https://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
% 
% -------------------------------------------------------------------------
close all
clear all
clc

% STL File Reader routines path:
%addpath('/YOUR_STLREAD_PATH/STLRead')

% Directory containing .stl file:
basedir='./';
% Name of .stl file:
file='dome.stl';

%% Parameters:
IAXIS = 1; JAXIS = 2; KAXIS = 3; MDIM = 3;
NOD1  = 1; NOD2  = 2; NOD3  = 3; MNOD = 3;
fast = true; % Use matlab functions in conversion stl -> FEM_MESH.
plot_normals    = false;
do_sanity_tests = false;

%% Import an STL mesh:
fprintf('%s \b','1. Reading binary .stl file...  ')
tstart=cputime;
fv      = stlread([basedir file]);
[F,V,N] = stlread([basedir file]);
nfaces=length(N(:,IAXIS));
telapsed=cputime-tstart;
fprintf('%s \n',['STL file read. ' num2str(telapsed) ' sec'])

%% Render Mesh:
fprintf('%s \b','2. Plotting mesh...  ')
tstart=cputime;
patch(fv,'FaceColor',       [0.8 0.8 1.0],        ...
         'FaceAlpha',       [0.5],                ...
         'EdgeColor',       [0. 0. 0.], ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
camlight('headlight');
material('dull');
axis('image');
view([-135 35]);

% Add normals to check on display that they are pointing correctly outside
% of the object.
if(plot_normals)
    hold on
    a=0.5;
    xyz1=1/3*(V(F(1:nfaces,NOD1),IAXIS:KAXIS) + ...
          V(F(1:nfaces,NOD2),IAXIS:KAXIS) + ...
          V(F(1:nfaces,NOD3),IAXIS:KAXIS));
    xyz2=xyz1(1:nfaces,IAXIS:KAXIS)+a*N(1:nfaces,IAXIS:KAXIS);     
    plot3([xyz1(1:nfaces,IAXIS) xyz2(1:nfaces,IAXIS)]', ...
          [xyz1(1:nfaces,JAXIS) xyz2(1:nfaces,JAXIS)]', ...
          [xyz1(1:nfaces,KAXIS) xyz2(1:nfaces,KAXIS)]','r');
end  
telapsed=cputime-tstart;
fprintf('%s \n',['Plot completed. ' num2str(telapsed) ' sec'])

%% Reduce STL format to FEM_MESH:
fprintf('%s \b','3. Converting STL format to FEM_MESH...  ')
tstart=cputime;
% Geometry threshold to identify different locations in space.
GEOMEPS = 10^-14*(max(max(V))-min(min(V)));
% Number of nodes obtained from .stl file:
stlnodes=length(V(:,IAXIS));
WSELEM=zeros(nfaces,MNOD);

if (~fast)
    IND=zeros(1,stlnodes);
    XYZ=zeros(stlnodes,MDIM);
    % First Node XYZ:
    nnodes=1;
    inod=1;
    XYZ(nnodes,IAXIS:KAXIS)=V(inod,IAXIS:KAXIS);
    IND(inod)=nnodes;
    % Fill rest of XYZ and IND:
    for inod=2:stlnodes
       inlist=0;
       for jnod=1:nnodes
           if(norm(V(inod,IAXIS:KAXIS)-XYZ(jnod,IAXIS:KAXIS)) < GEOMEPS) 
               inlist=1;
               IND(inod)=jnod;
               break
           end
       end
       if(~inlist)
           nnodes=nnodes+1;
           XYZ(nnodes,IAXIS:KAXIS)=V(inod,IAXIS:KAXIS);
           IND(inod)=nnodes;
       end
    end
    XYZ=XYZ(1:nnodes,IAXIS:KAXIS);

    % Re-index faces:
    for iface=1:nfaces
        WSELEM(iface,NOD1:NOD3)=IND(F(iface,NOD1:NOD3));
    end

else
  
    [XYZ,IA,IC]=uniquetol(V,GEOMEPS,'ByRows',true,'DataScale',1);
    nnodes = length(XYZ(:,IAXIS));
    % Re-index faces:
    for iface=1:nfaces
        WSELEM(iface,NOD1:NOD3)=IC(F(iface,NOD1:NOD3));
    end
    
end
telapsed=cputime-tstart;
fprintf('%s \n',['Conversion finished. ' num2str(telapsed) ' sec'])

%% Sanity Tests:
fprintf('%s \b','4. Test for repeated nodes, faces...  ')
tstart=cputime;
if(do_sanity_tests)
    % Check we don't have repeated nodes:
    for inod=1:nnodes
        DIFF=abs(XYZ(:,IAXIS)-XYZ(inod,IAXIS)) + ...
             abs(XYZ(:,JAXIS)-XYZ(inod,JAXIS)) + ...
             abs(XYZ(:,KAXIS)-XYZ(inod,KAXIS));
        [vdf,idf]=find(DIFF < GEOMEPS);
        if (length(vdf)>1)
            disp(['Found same nodes: ' num2str([inod idf])])
        end
    end

    % Check we don't have repeated faces:
    for iface=1:nfaces
        DIFF=abs(WSELEM(:,NOD1)-WSELEM(iface,NOD1)) + ...
             abs(WSELEM(:,NOD2)-WSELEM(iface,NOD2)) + ...
             abs(WSELEM(:,NOD3)-WSELEM(iface,NOD3));
        [vdf,idf]=find(DIFF < GEOMEPS);
        if (length(vdf)>1)
            disp(['Found same faces: ' num2str([iface idf])])
        end
    end
end

telapsed=cputime-tstart;
if(do_sanity_tests)
    fprintf('%s \n\n',['Tests finished. ' num2str(telapsed) ' sec'])
else
    fprintf('%s \n\n',['Tests skipped. ' num2str(telapsed) ' sec'])
end
fprintf('%s \n\n',['FEM_MESH contains ' num2str(nnodes) ' VERTS and ' ...
                                        num2str(nfaces) ' FACES.'])

%% Finally write geometry file with same name as file and extension .dat.
% The contents of this ASCII file are to be copied to the FDS input file.
fprintf('%s \b','5. Write FEM_MESH geometry in ASCII format...  ')
tstart=cputime;
fileout=[file(1:end-3) 'dat'];
[fid]=fopen([basedir fileout],'w');

GEOM_ID='FEM_MESH';
SURF_ID=file(1:end-4);

% Write &GEOM namelist
[wid]=fprintf(fid,'&GEOM ID=''%8s'', SURF_ID=''%5s''\n',GEOM_ID,SURF_ID);

% Vertices:
[wid]=fprintf(fid,'VERTS=\n');
for inod=1:nnodes
    [wid]=fprintf(fid,' %14.8f,   %14.8f,   %14.8f,\n',XYZ(inod,IAXIS:KAXIS));
end

% Tetrahedra:
% No volume elements.

% Surface Triangles:
[wid]=fprintf(fid,'FACES=\n');
for iwsel=1:nfaces
    [wid]=fprintf(fid,' %6d,   %6d,   %6d, %6d\n',WSELEM(iwsel,NOD1:NOD3),1);
end

[wid]=fprintf(fid,'/ \n');

% Close geom file
fclose(fid);
telapsed=cputime-tstart;
fprintf('%s \n\n',['File written. ' num2str(telapsed) ' sec'])
fprintf('%s \n\n',['File ' basedir fileout ' written. '])

return
