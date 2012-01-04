% McDermott
% 9-19-11
% sphere_demo.m

close all
clear all

n=32;
r = 1;
theta = linspace(0,pi,n);
phi = linspace(0,2*pi,n);

x = zeros(n^2,1);
y = zeros(n^2,1);
z = zeros(n^2,1);

p=0;
for i=1:n
    for j=1:n
        p = p+1;
        x(p) = r*sin(theta(i))*cos(phi(j));
        y(p) = r*sin(theta(i))*sin(phi(j));
        z(p) = r*cos(theta(i));
    end
end

F = convhulln([x,y,z]);
trisurf(F,x,y,z)

% write FDS input file

fid = fopen('sphere_demo.fds','wt');

head = ['&HEAD CHID=''sphere_demo'', TITLE=''Demo of Matlab convhulln, created by sphere_demo.m'' /'];
fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['&MESH IJK=64,32,32, XB=-4,4,-2,2,-2,2 /']; fprintf(fid,'%s\n',mesh);

fprintf(fid,'%s\n','  '); % blank line

time = ['&TIME T_END=10./']; fprintf(fid,'%s\n',time);

fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=0 /']; fprintf(fid,'%s\n',misc);

fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT MB=''XMIN'', SURF_ID=''supply'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''supply'', VEL=-1., COLOR=''LIGHT BLUE''/'];  fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''sphere'', COLOR=''GREEN''/'];  fprintf(fid,'%s\n',surf);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBY=0, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-4,4,-2,2,-2,2, QUANTITY=''P MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-4,4,-2,2,-2,2, QUANTITY=''U MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-4,4,-2,2,-2,2, QUANTITY=''V MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-4,4,-2,2,-2,2, QUANTITY=''W MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

% write VERT lines

for i=1:length(x)
    vert = ['&VERT X=',num2str(x(i)),',',num2str(y(i)),',',num2str(z(i)),' /']; fprintf(fid,'%s\n',vert);
end

fprintf(fid,'%s\n','  '); % blank line

% write FACE lines

surf_id='''sphere''';
for i=1:length(F(:,1))
    face = ['&FACE N=',num2str(F(i,1)),',',num2str(F(i,3)),',',num2str(F(i,2)),', SURF_ID=',surf_id,' /']; fprintf(fid,'%s\n',face);
end

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    