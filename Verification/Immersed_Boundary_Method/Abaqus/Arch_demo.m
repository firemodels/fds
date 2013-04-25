
close all
clear all

% read the Abaqus input file

[X Tet] = abaqus2fds('Arch_demo.inp');

% write FDS input file

fid = fopen('Arch_demo2.fds','wt');

head = ['&HEAD CHID=''Arch_demo'', TITLE=''Demo of Abaqus >> Matlab freeBoundary, created by Arch_demo.m'' /']; fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['&MESH IJK=32,32,32, XB=-1.,1.,-1.,1.,0,1/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=60./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=1 /']; fprintf(fid,'%s\n',misc); fprintf(fid,'%s\n','  ');

dump = ['&DUMP CUTCELL_DATA_FILE=.TRUE. /']; fprintf(fid,'%s\n',dump);

fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''hot vent'', COLOR=''RED'', CONVECTIVE_HEAT_FLUX=10./'];  fprintf(fid,'%s\n',surf);
part = ['&PART ID=''p1'', SURF_ID=''arch surf'', STATIC=.TRUE., PROP_ID=''ball'', QUANTITIES(1)=''PARTICLE TEMPERATURE'' /']; fprintf(fid,'%s\n',part);
prop = ['&PROP ID=''ball'', SMOKEVIEW_ID=''SPHERE'', SMOKEVIEW_PARAMETERS(1)=''D=0.01'' /']; fprintf(fid,'%s\n',prop);
surf = ['&SURF ID=''arch surf'', COLOR=''BLUE'', EMISSIVITY=0.95, RADIUS=0.005, GEOMETRY=''SPHERICAL'' /',];  fprintf(fid,'%s\n',surf);

fprintf(fid,'%s\n','  '); % blank line

geom = ['&GEOM ID=''arch'', SHAPE=''COMPLEX'', SURF_ID=''cold'',DT_GEOC=1., GEOC_FILENAME=''Arch_demo_bending_01.gc'' /']; fprintf(fid,'%s\n',geom);

fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT XB=-.2,.2,-.2,.2,0,0, SURF_ID=''hot vent''/'];  fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

matl = ['&MATL ID=''my solid'', DENSITY=2700., CONDUCTIVITY=120., SPECIFIC_HEAT=.903 /'];  fprintf(fid,'%s\n',matl);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBY=0, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=0, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

fprintf(fid,'%s\n','  '); % blank line

bnde = ['&BNDE QUANTITY=''WALL TEMPERATURE'' /'];  fprintf(fid,'%s\n',bnde);  fprintf(fid,'%s\n','  '); % blank line

% write VERT lines

for i=1:length(X(:,1))
    vert = ['&VERT X=',num2str(X(i,1)),',',num2str(X(i,2)),',',num2str(X(i,3)),' /']; fprintf(fid,'%s\n',vert);
end

fprintf(fid,'%s\n','  '); % blank line

% write FACE lines

TR = TriRep(Tet,X);
FF = freeBoundary(TR); % freeBoundary uses the same ordering convention as FDS: right-hand rule for outward normal

for i=1:length(FF(:,1))
    surf_id='''arch surf''';
    part_id='''p1''';
    face = ['&FACE N=',num2str(FF(i,1)),',',num2str(FF(i,2)),',',num2str(FF(i,3)),', SURF_ID=',surf_id, ',PART_ID=', part_id, ' /' ]; fprintf(fid,'%s\n',face);
end

fprintf(fid,'%s\n','  '); % blank line

% % write VOLU lines
% 
% for i=1:length(Tet(:,1))
%     volu = ['&VOLU N=',num2str(Tet(i,1)),',',num2str(Tet(i,2)),',',num2str(Tet(i,3)),',',num2str(Tet(i,4)),',', ...
%         ' MATL_ID=''my solid'' /'];
%     
%     fprintf(fid,'%s\n',volu);
% end
% 
% fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);



    