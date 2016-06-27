% McDermott
% 5-30-14
% sphere_fire_demo.m

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
        z(p) = r*cos(theta(i)) + 2;
    end
end

F = convhulln([x,y,z]);
trisurf(F,x,y,z)

% write FDS input file

fid = fopen('sphere_fire_demo.fds','wt');

head = ['&HEAD CHID=''sphere_fire_demo'', TITLE=''Demo of Matlab convhulln, created by sphere_fire_demo.m'' /'];
fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['&MESH IJK=32,32,32, XB=-2,2,-2,2,0,4 /']; fprintf(fid,'%s\n',mesh);

fprintf(fid,'%s\n','  '); % blank line

time = ['&TIME T_END=10./']; fprintf(fid,'%s\n',time);

fprintf(fid,'%s\n','  '); % blank line

reac = ['&REAC FUEL=''PROPANE'', SOOT_YIELD=0.01/'];  fprintf(fid,'%s\n',reac);

fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT XB=-1,1,-1,1,0,0, SURF_ID=''fire'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''fire'', HRRPUA=1000., COLOR=''RED''/'];  fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''sphere'', COLOR=''GREEN''/'];  fprintf(fid,'%s\n',surf);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBY=0, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=0, QUANTITY=''TEMPERATURE'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=0, QUANTITY=''DIVERGENCE'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=0, QUANTITY=''MASS FRACTION'', SPEC_ID=''PROPANE'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

bnde = ['&BNDE QUANTITY=''GAS TEMPERATURE'' /'];  fprintf(fid,'%s\n',bnde);

fprintf(fid,'%s\n','  '); % blank line

% write GEOM lines

nx = length(x);
geom  = ['&GEOM ID=''sphere'', SURF_ID=''sphere'',']; fprintf(fid,'%s\n',geom);
verts = ['      VERTS=',num2str(x(1)),',',num2str(y(1)),',',num2str(z(1)),',']; fprintf(fid,'%s\n',verts);
for i=2:nx-1
verts = ['            ',num2str(x(i)),',',num2str(y(i)),',',num2str(z(i)),',']; fprintf(fid,'%s\n',verts);
end
verts = ['            ',num2str(x(nx)),',',num2str(y(nx)),',',num2str(z(nx)),',']; fprintf(fid,'%s\n',verts);

% write FACES

nf = length(F(:,1));
faces = ['      FACES=',num2str(F(1,1)),',',num2str(F(1,3)),',',num2str(F(1,2)),',']; fprintf(fid,'%s\n',faces);
for i=2:nf-1
faces = ['            ',num2str(F(i,1)),',',num2str(F(i,3)),',',num2str(F(i,2)),',']; fprintf(fid,'%s\n',faces);
end
faces = ['            ',num2str(F(nf,1)),',',num2str(F(nf,3)),',',num2str(F(nf,2)),'/']; fprintf(fid,'%s\n',faces);

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    