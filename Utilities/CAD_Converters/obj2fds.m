%% ------------------------------------------------------------------------
% Script obj2fds:
%
% This script reads a Wavefront .obj CAD file and writes the geometry nodes 
% and surface triangles in FDS GEOM format. The result can be
% directly copied to FDS input files.
%
% In order to use this script you need to download the script read_wobj.m
% -------------------------------------------------------------------------
close all
clear all
clc

file='Death_Star.obj';
%file='UWO-BLWT-SS21-Test6-40ft.obj';

%% Parameters:
IAXIS = 1; JAXIS = 2; KAXIS = 3; MDIM = 3;
NOD1  = 1; NOD2  = 2; NOD3  = 3; MNOD = 3;

OBJ = read_wobj(file);

V = OBJ.vertices;
F = OBJ.objects(3).data.vertices;
N = OBJ.vertices_normal;

nnodes = length(V(:,IAXIS));
nfaces = 0;
nsurf = 0;
SURF_INDEX = zeros([1 20]);
SURF_ID = strings([1 20]);

% Draw the object

FV.vertices=OBJ.vertices;
figure
for i=1:length(OBJ.objects)
    if strcmp(OBJ.objects(i).type,'usemtl')
        flag = 0;
        for ii=1:i-1
            if (strcmp(SURF_ID(ii),OBJ.objects(i).data)); flag=1; end
        end      
        if flag==0
            nsurf = nsurf + 1;
            SURF_ID(nsurf) = OBJ.objects(i).data;
            SURF_INDEX(i+1) = nsurf;
        end
    end
    if strcmp(OBJ.objects(i).type,'f')
        FV.faces=OBJ.objects(i).data.vertices;
        nfaces = nfaces + length(FV.faces(:,NOD1));
        patch(FV,'facecolor',[1 0 0]); camlight ;hold on
    end
end

fprintf('%s \b','Converting OBJ format to FDS GEOM...  ')

fprintf('%s \n\n',['FEM_MESH contains ' num2str(nnodes) ' VERTS and ' ...
                                        num2str(nfaces) ' FACES.'])

fileout=[file(1:end-3) 'dat'];
[fid]=fopen([fileout],'w');

label{1} = ['&GEOM ID=''%8s'', ']; 
for i=2:10
    label{i} = ['SURF_ID(' num2str(i-1) ')=''%5s'', '];
end
formatSpec = [strcat(label{1:nsurf+1}) '/\n'];
GEOM_ID='FEM_MESH';

title_string = sprintf(formatSpec,GEOM_ID,SURF_ID(1:nsurf));
[wid]=fprintf(fid,title_string);

[wid]=fprintf(fid,'VERTS=\n');
for inod=1:nnodes
    [wid]=fprintf(fid,' %14.8f,   %14.8f,   %14.8f,\n',V(inod,IAXIS:KAXIS));
end

[wid]=fprintf(fid,'FACES=\n');
for i=1:length(OBJ.objects)
    if strcmp(OBJ.objects(i).type,'f')
        FV.faces=OBJ.objects(i).data.vertices;
        for ii=1:length(FV.faces(:,NOD1))
            [wid]=fprintf(fid,' %6d,%6d,%6d,%6d\n',FV.faces(ii,NOD1:NOD3),SURF_INDEX(i));
        end
    end
end

[wid]=fprintf(fid,'/ \n');

fclose(fid);

fprintf('%s \n\n',['File ' fileout ' written. '])

return
