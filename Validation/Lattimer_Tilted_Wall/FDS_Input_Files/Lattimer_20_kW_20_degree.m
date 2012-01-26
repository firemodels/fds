% McDermott
% 10-19-11
% Arch_demo.m
%
% Builds I-beam geometry from tet volume mesh created in Abaqus.
% The surface tri mesh is generated using the freeBoundary function.

close all
clear all

% read the Abaqus input file

[X Tet] = importabaqus('Lattimer_20_kW_20_degree.inp');

% write FDS input file

fid = fopen('Lattimer_20_kW_20_degree.fds','wt');

head = ['&HEAD CHID=''Lattimer_20_kW_20_degree'', TITLE=''Surface temperature on slopped wall, 20 kW propane fire'' /'];
fprintf(fid,'%s\n',head);

fprintf(fid,'%s\n','  '); % blank line

mesh = ['MESH IJK=30,30,30, XB=-0.75,0.75,-0.20,1.30,-0.45,1.05 /']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');
mesh = ['&MULT ID=''m1'', DX=.5, DY=.5, DZ=.5, I_LOWER=-1,I_UPPER=1, J_UPPER=2, K_LOWER=-1,K_UPPER=1/']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');
mesh = ['&MESH IJK=10,10,10, XB=-0.25,0.25,-0.20,0.30,0.05,0.55 MULT_ID=''m1'' /']; fprintf(fid,'%s\n',mesh); fprintf(fid,'%s\n','  ');

time = ['&TIME T_END=1000./']; fprintf(fid,'%s\n',time); fprintf(fid,'%s\n','  '); % blank line

dump = ['&DUMP DT_SLCF=0.1, DT_BNDF=1., DT_DEVC_LINE=100., DT_DEVC=1. /']; fprintf(fid,'%s\n',dump);

fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC IMMERSED_BOUNDARY_METHOD=1 /']; fprintf(fid,'%s\n',misc);
misc = ['&MISC BNDF_DEFAULT=.FALSE. /']; fprintf(fid,'%s\n',misc);
misc = ['&MISC GVEC=0.,-9.81,0. /']; fprintf(fid,'%s\n',misc);

fprintf(fid,'%s\n','  '); % blank line

devc = ['&DEVC XYZ=0,0,0.3, ID=''burner on'', QUANTITY=''TIME'', SETPOINT=900., INITIAL_STATE=.TRUE. /']; fprintf(fid,'%s\n',devc);
surf = ['&SURF ID=''burner'', HRRPUA=222., TMP_FRONT=400., COLOR=''RED'' /']; fprintf(fid,'%s\n',surf);
vent = ['&VENT XB=0.00,0.30,0.001,0.001,0.15,0.45, SURF_ID=''burner'', DEVC_ID=''burner on'' /']; fprintf(fid,'%s\n',vent);
obst = ['&OBST XB=0.00,0.30,-0.1,0.00,0.15,0.45 /']; fprintf(fid,'%s\n',obst);

fprintf(fid,'%s\n','  '); % blank line

matl = ['&MATL ID=''ceramic'', CONDUCTIVITY=.15, SPECIFIC_HEAT=1., DENSITY=272. /']; fprintf(fid,'%s\n',matl);
surf = ['&SURF ID=''ceramic board'', EMISSIVITY=0.95, COLOR=''BLUE'', MATL_ID=''ceramic'', THICKNESS=0.0254, BACKING=''EXPOSED'' /']; fprintf(fid,'%s\n',surf);
%obst = ['&OBST XB=-0.0254,-0.0,0,0.915,-0.005,0.605, SURF_ID=''ceramic board'', BNDF_FACE(1)=.TRUE. /']; fprintf(fid,'%s\n',obst);

fprintf(fid,'%s\n','  '); % blank line

%surf = ['&SURF ID=''hot'', COLOR=''RED'', CONVECTIVE_HEAT_FLUX=1000./'];  fprintf(fid,'%s\n',surf);
reac = ['&REAC FUEL=''PROPANE'', SOOT_YIELD=0.01/'];  fprintf(fid,'%s\n',reac);

fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT MB=''XMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMIN'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'', SURF_ID=''OPEN''/'];  fprintf(fid,'%s\n',vent);

fprintf(fid,'%s\n','  '); % blank line

bndf = ['&BNDF QUANTITY=''GAS TEMPERATURE'' /'];  fprintf(fid,'%s\n',bndf);

fprintf(fid,'%s\n','  '); % blank line

% devc = ['&DEVC XB=0,0,0,0.8,0.3,0.3, QUANTITY=''WALL TEMPERATURE'', ID=''TIR'', IOR=1, POINTS=40, DEVC_ID=''burner on'' /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.075,0.3, QUANTITY=''WALL TEMPERATURE'', ID=''TC1'', IOR=1 /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.299,0.3, QUANTITY=''WALL TEMPERATURE'', ID=''TC2'', IOR=1 /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.600,0.3, QUANTITY=''WALL TEMPERATURE'', ID=''TC3'', IOR=1 /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.075,0.3, QUANTITY=''ADIABATIC SURFACE TEMPERATURE'', ID=''AST1'', IOR=1 /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.075,0.3, QUANTITY=''ADIABATIC SURFACE TEMPERATURE'', ID=''AST2'', IOR=1 /']; fprintf(fid,'%s\n',devc);
% devc = ['&DEVC XYZ=0,0.075,0.3, QUANTITY=''ADIABATIC SURFACE TEMPERATURE'', ID=''AST3'', IOR=1 /']; fprintf(fid,'%s\n',devc);

fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF PBZ=0.3, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBZ=0.3, QUANTITY=''TEMPERATURE'', VECTOR=.TRUE. /'];  fprintf(fid,'%s\n',slcf);

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
    surf_id='''ceramic board''';
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



    