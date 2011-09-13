% McDermott
% 9-9-11
% peaks_demo.m

close all
clear all

[x,y] = meshgrid(1:15,1:15);
tri = delaunay(x,y);
z = peaks(15);
%trimesh(tri,x,y,z)

% write FDS input file

fid = fopen('peaks_demo.fds','wt');

head = ['&HEAD CHID=''peaks_demo'', TITLE=''Matlab demo of delaunay'' /']; fprintf(fid,'%s\n',head); fprintf(fid,'%s\n','  ');

mesh = ['&MESH IJK=30,30,30, XB=0,15,0,15,-7.5,7.5/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=0./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''terrain'', COLOR=''GREEN''/']; fprintf(fid,'%s\n',surf); fprintf(fid,'%s\n','  ');

vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMIN'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN''/']; fprintf(fid,'%s\n',vent);

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



    