% McDermott
% 10-19-11
% Ibeam_demo.m
%
% Builds I-beam geometry from tet volume mesh created in Abaqus.
% The surface tri mesh is generated using the freeBoundary function.

close all
clear all

% read the Abaqus input file

[X Tet] = abaqus2fds('Ibeam.inp');

% write FDS input file

fid = fopen('Ibeam_demo.fds','wt');

head = ['&HEAD CHID=''Ibeam_demo'', TITLE=''Demo of Abaqus >> Matlab freeBoundary, created by Ibeam_demo.m'' /']; fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['&MESH IJK=32,32,32, XB=-.5,.5,-.5,.5,0,1/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=10./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=1 /']; fprintf(fid,'%s\n',misc);

fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''hot'', COLOR=''RED'', CONVECTIVE_HEAT_FLUX=1000./'];  fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''cold'', COLOR=''BLUE''/'];  fprintf(fid,'%s\n',surf);

fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);

vent = ['&VENT XB=-.2,.2,-.2,.2,0,0, SURF_ID=''hot''/'];  fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

matl = ['&MATL ID=''my solid'', DENSITY=1000., CONDUCTIVITY=1., SPECIFIC_HEAT=1. /'];  fprintf(fid,'%s\n',matl);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBY=0, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=0, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

bnde = ['&BNDE QUANTITY=''WALL TEMPERATURE'' /'];  fprintf(fid,'%s\n',bnde);

fprintf(fid,'%s\n','  '); % blank line

% write VERT lines

for i=1:length(X(:,1))
    vert = ['&VERT X=',num2str(X(i,1)),',',num2str(X(i,2)),',',num2str(X(i,3)),' /']; fprintf(fid,'%s\n',vert);
end

fprintf(fid,'%s\n','  '); % blank line

% write FACE lines

TR = TriRep(Tet,X);
FF = freeBoundary(TR); % freeBoundary uses the same ordering convention as FDS: right-hand rule for outward normal

for i=1:length(FF(:,1))
    surf_id='''cold''';
    face = ['&FACE N=',num2str(FF(i,1)),',',num2str(FF(i,2)),',',num2str(FF(i,3)),', SURF_ID=',surf_id,' /']; fprintf(fid,'%s\n',face);
end

fprintf(fid,'%s\n','  '); % blank line

% write VOLU lines

for i=1:length(Tet(:,1))
    volu = ['&VOLU N=',num2str(Tet(i,1)),',',num2str(Tet(i,2)),',',num2str(Tet(i,3)),',',num2str(Tet(i,4)),',', ...
        ' MATL_ID=''my solid'' /'];
    
    fprintf(fid,'%s\n',volu);
end

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    