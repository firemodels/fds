% McDermott
% 9-13-11
% cube_demo.m
%
% Builds cube with tet volume mesh in FDS format from Matlab tetramesh example.
% The surface tri mesh is generated using the freeBoundary function.

close all
clear all

d = [-1 1];
[x,y,z] = meshgrid(d,d,d);  % A cube
x = [x(:);0];
y = [y(:);0];
z = [z(:);0];    % [x,y,z] are corners of a cube plus the center.
dt = DelaunayTri(x,y,z);
Tes = dt(:,:);
X = [x(:) y(:) z(:)];
%tetramesh(Tes,X);camorbit(20,0)

% write FDS input file

fid = fopen('cube_demo.fds','wt');

head = ['&HEAD CHID=''cube_demo'', TITLE=''Demo of Matlab tetramesh and freeBoundary, created by cube_demo.m'' /'];
fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['&MESH IJK=16,16,16, XB=-5,5,-5,5,-5,5/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=10./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=0 /']; fprintf(fid,'%s\n',misc); fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT XB=-5,5,-5,5,-5,-5, SURF_ID=''supply'' /']; fprintf(fid,'%s\n',vent);
vent = ['&VENT XB=-5,5,-5,5,5,5, SURF_ID=''OPEN'' /']; fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''supply'', VEL=-1., COLOR=''LIGHT BLUE''/'];  fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s1'', COLOR=''WHITE''/'];  fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s2'', COLOR=''YELLOW''/']; fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s3'', COLOR=''ORANGE''/']; fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s4'', COLOR=''RED''/'];    fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s5'', COLOR=''BLUE''/'];   fprintf(fid,'%s\n',surf);
surf = ['&SURF ID=''s6'', COLOR=''GREEN''/'];  fprintf(fid,'%s\n',surf);

fprintf(fid,'%s\n','  '); % blank line

matl = ['&MATL ID=''my solid'', DENSITY=1000., CONDUCTIVITY=1., SPECIFIC_HEAT=1. /'];  fprintf(fid,'%s\n',matl);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBY=0, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-5,5,-5,5,-5,5, QUANTITY=''P MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-5,5,-5,5,-5,5, QUANTITY=''U MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-5,5,-5,5,-5,5, QUANTITY=''V MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF XB=-5,5,-5,5,-5,5, QUANTITY=''W MASK'', CELL_CENTERED=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

% write VERT lines

for i=1:length(x)
    vert = ['&VERT X=',num2str(x(i)),',',num2str(y(i)),',',num2str(z(i)),' /']; fprintf(fid,'%s\n',vert);
end

fprintf(fid,'%s\n','  '); % blank line

% write FACE lines

TR = TriRep(Tes,x,y,z);
FF = freeBoundary(TR); % freeBoundary uses the same ordering convention as FDS: right-hand rule for outward normal

for i=1:length(FF(:,1))
    
    % compute face normal (only needed here so we can assign SURF_ID based on orientation)
    u(1) = x(FF(i,3))-x(FF(i,1));
    u(2) = y(FF(i,3))-y(FF(i,1));
    u(3) = z(FF(i,3))-z(FF(i,1));
    v(1) = x(FF(i,2))-x(FF(i,1));
    v(2) = y(FF(i,2))-y(FF(i,1));
    v(3) = z(FF(i,2))-z(FF(i,1));
    
    FN(1) = u(2)*v(3)-u(3)*v(2);
    FN(2) = u(3)*v(1)-u(1)*v(3);
    FN(3) = u(1)*v(2)-u(2)*v(1);
    
    LENGTH = sqrt(sum(FN.*FN));
    FN = FN/LENGTH;
    
    if FN(1)==-1; surf_id='''s1'''; end
    if FN(1)==1;  surf_id='''s2'''; end
    
    if FN(2)==-1; surf_id='''s3'''; end
    if FN(2)==1;  surf_id='''s4'''; end
    
    if FN(3)==-1; surf_id='''s5'''; end
    if FN(3)==1;  surf_id='''s6'''; end
    
    face = ['&FACE N=',num2str(FF(i,1)),',',num2str(FF(i,2)),',',num2str(FF(i,3)),', SURF_ID=',surf_id,' /']; fprintf(fid,'%s\n',face);
end

fprintf(fid,'%s\n','  '); % blank line

% write VOLU lines

for i=1:length(dt(:,1))
    volu = ['&VOLU N=',num2str(dt.Triangulation(i,1)),',', ...
        num2str(dt.Triangulation(i,2)),',', ...
        num2str(dt.Triangulation(i,3)),',', ...
        num2str(dt.Triangulation(i,4)),',', ...
        ' MATL_ID=''my solid'' /'];
    
    fprintf(fid,'%s\n',volu);
end

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    