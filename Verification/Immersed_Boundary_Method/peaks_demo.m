% McDermott
% 9-9-11
% peaks_demo.m

close all
clear all

N=16;
[x,y] = meshgrid(0:N,0:N);
tri = delaunay(x,y);
z = peaks(N+1);
%trimesh(tri,x,y,z)

% write FDS input file

fid = fopen('peaks_demo.fds','wt');

head = ['&HEAD CHID=''peaks_demo'', TITLE=''Matlab demo of delaunay'' /']; fprintf(fid,'%s\n',head); fprintf(fid,'%s\n','  ');

mesh = ['&MESH IJK=32,32,32, XB=0,16,0,16,-8,8/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=1./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=0 /']; fprintf(fid,'%s\n',misc); fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''terrain'', COLOR=''GREEN''/']; fprintf(fid,'%s\n',surf); fprintf(fid,'%s\n','  ');

vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF XB=0,16,0,16,-8,8, QUANTITY=''P MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=0,16,0,16,-8,8, QUANTITY=''U MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=0,16,0,16,-8,8, QUANTITY=''V MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=0,16,0,16,-8,8, QUANTITY=''W MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

% write VERT lines

s = size(x);

for i=1:s(1)
    for j=1:s(2)
       vert = ['&VERT X=',num2str(x(i,j)),',',num2str(y(i,j)),',',num2str(z(i,j)),' /'];
       fprintf(fid,'%s\n',vert);
    end
end

fprintf(fid,'%s\n','  '); % blank line

% write FACE lines

for i=1:length(tri(:,1))
    face = ['&FACE N=',num2str(tri(i,1)),',',num2str(tri(i,3)),',',num2str(tri(i,2)),', SURF_ID=''terrain'' /'];
    fprintf(fid,'%s\n',face);
end

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    