% McDermott
% 8-31-09
% Askervein_topo.m

close all
clear all
fclose('all');

% terrain data

x_data = csvread('x_data.csv');
y_data = csvread('y_data.csv');
Z_data = csvread('z_data.csv');

x_data = x_data - min(x_data);
y_data = y_data - min(y_data);

surf(x_data,y_data,Z_data)
xlabel('x')
ylabel('y')
zlabel('z')

min_x_data = min(x_data);
max_x_data = max(x_data);
min_y_data = min(y_data);
max_y_data = max(y_data);
x_min = min_x_data - 0.05*(max_x_data-min_x_data);
x_max = max_x_data + 0.05*(max_x_data-min_x_data);
y_min = min_y_data - 0.05*(max_y_data-min_y_data);
y_max = max_y_data + 0.05*(max_y_data-min_y_data);

min_z_data = min(min(Z_data));
z_min = min_z_data - 20;
z_max = 500; % based on Adam's presentation

% set FDS input file parameters

T = 3600;
% N = [128 128 64]
% dx = (x_max-x_min)/N(1)
% dy = (y_max-y_min)/N(2)
% dz = (z_max-z_min)/N(3)
dx = 20;
dy = 20;
dz = 10;
N = [round((x_max-x_min)/dx), round((y_max-y_min)/dy), round((z_max-z_min)/dz)]
for i=1:3
    if mod(N(i),2)~=0; N(i)=N(i)+1; end
end
x_max = x_min+N(1)*dx;
y_max = y_min+N(2)*dy;
z_max = z_min+N(3)*dz;

xf = x_min:dx:x_min+N(1)*dx;
yf = y_min:dy:y_min+N(2)*dy;

[X,Y] = meshgrid(x_data,y_data);

% create the FDS input file

fid = fopen('Askervein.fds','wt');

head = ['&HEAD CHID=''Askervein'', TITLE=''Wind flow over Askervein hill'' /'];
fprintf(fid,'%s\n',head);
fprintf(fid,'%s\n','  '); % blank line

% place fine mesh first for DEVC output
n1 = num2str(N(1));
n2 = num2str(N(2));
n3 = num2str(N(3));
x1 = num2str(x_min,'%1.2f');
x2 = num2str(x_max,'%1.2f');
y1 = num2str(y_min,'%1.2f');
y2 = num2str(y_max,'%1.2f');
z1 = num2str(z_min,'%1.2f');
z2 = num2str(z_max,'%1.2f');
mesh = ['&MESH IJK=',n1,',',n2,',',n3,', XB=',x1,',',x2,',',y1,',',y2,',',z1,',',z2,'/']; fprintf(fid,'%s\n',mesh);
% dx2 = 2*dx;
% dy2 = 2*dy;
% dz2 = 2*dz;
% N2(1) = (x_max-x_min)/dx2 + 2*floor(0.5*(x_max-x_min)/dx2);
% N2(2) = (y_max-y_min)/dy2 + 2*floor(0.5*(y_max-y_min)/dy2);
% N2(3) = (z_max-z_min)/dz2;
% n1 = num2str(N2(1));
% n2 = num2str(N2(2));
% n3 = num2str(N2(3));
% x1 = num2str(x_min-floor(0.5*(x_max-x_min)/dx2)*dx2,'%1.2f');
% x2 = num2str(x_max+floor(0.5*(x_max-x_min)/dx2)*dx2,'%1.2f');
% y1 = num2str(y_min-floor(0.5*(y_max-y_min)/dy2)*dy2,'%1.2f');
% y2 = num2str(y_max+floor(0.5*(y_max-y_min)/dy2)*dy2,'%1.2f');
% z1 = num2str(z_min,'%1.2f');
% z2 = num2str(z_max,'%1.2f');
% mesh = ['&MESH IJK=',n1,',',n2,',',n3,', XB=',x1,',',x2,',',y1,',',y2,',',z1,',',z2,'/']; fprintf(fid,'%s\n',mesh);
fprintf(fid,'%s\n','  '); % blank line

%trnz = ['&TRNZ IDERIV=0, CC=250, PC=100, MESH_NUMBER=1/']; fprintf(fid,'%s\n',trnz);
%trnz = ['&TRNZ IDERIV=1, CC=250, PC=.8, MESH_NUMBER=1/']; fprintf(fid,'%s\n',trnz);
%trnz = ['&TRNZ IDERIV=0, CC=250, PC=100, MESH_NUMBER=2/']; fprintf(fid,'%s\n',trnz);
%trnz = ['&TRNZ IDERIV=1, CC=250, PC=.8, MESH_NUMBER=2/']; fprintf(fid,'%s\n',trnz);
% fprintf(fid,'%s\n','  '); % blank line

time = ['&TIME T_END=',num2str(T,'%1.1f'),'/'];
fprintf(fid,'%s\n',time);
fprintf(fid,'%s\n','  '); % blank line

dump = ['&DUMP NFRAMES=3600/'];
fprintf(fid,'%s\n',dump);
fprintf(fid,'%s\n','  '); % blank line

misc = ['&MISC TMPA=15.0'];         fprintf(fid,'%s\n',misc);
misc = ['      LAPSE_RATE=-0.01'];  fprintf(fid,'%s\n',misc);
misc = ['      HUMIDITY=95.0'];     fprintf(fid,'%s\n',misc);
% misc = ['      U0=6.0, FORCE_VECTOR(1)=0.01']; fprintf(fid,'%s\n',misc);
% misc = ['      V0=16.0,FORCE_VECTOR(2)=0.02/']; fprintf(fid,'%s\n',misc);
misc = ['      DT_MEAN_FORCING=60.0'];     fprintf(fid,'%s\n',misc);
misc = ['      U0=6.0, MEAN_FORCING(1)=T']; fprintf(fid,'%s\n',misc);
misc = ['      V0=16.0,MEAN_FORCING(2)=T/']; fprintf(fid,'%s\n',misc);
fprintf(fid,'%s\n','  '); % blank line

surf = ['&SURF ID=''terrain'',ROUGHNESS=0.02, DEFAULT=T/']; fprintf(fid,'%s\n',surf);
fprintf(fid,'%s\n','  '); % blank line

vent = ['&VENT MB=''XMIN'',SURF_ID=''OPEN'',WIND=T/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''XMAX'',SURF_ID=''OPEN'',WIND=T/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMIN'',SURF_ID=''OPEN'',WIND=T/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''YMAX'',SURF_ID=''OPEN'',WIND=T/']; fprintf(fid,'%s\n',vent);
vent = ['&VENT MB=''ZMAX'',SURF_ID=''OPEN'',WIND=T/']; fprintf(fid,'%s\n',vent);
fprintf(fid,'%s\n','  '); % blank line

slcf = ['&SLCF QUANTITY=''VELOCITY'',VECTOR=.TRUE.,AGL_SLICE=10.0/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=0,QUANTITY=''VELOCITY'',VECTOR=.TRUE./']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=0,QUANTITY=''VELOCITY'',VECTOR=.TRUE./']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=984,QUANTITY=''VELOCITY'',VECTOR=.TRUE./']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=1625,QUANTITY=''VELOCITY'',VECTOR=.TRUE./']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBZ=150,QUANTITY=''VELOCITY'',VECTOR=.TRUE./']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=984,QUANTITY=''TEMPERATURE''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=1625,QUANTITY=''TEMPERATURE''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBZ=150,QUANTITY=''TEMPERATURE''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=984,QUANTITY=''H''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=1625,QUANTITY=''H''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBZ=150,QUANTITY=''H''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBX=984,QUANTITY=''PRESSURE''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBY=1625,QUANTITY=''PRESSURE''/']; fprintf(fid,'%s\n',slcf);
slcf = ['&SLCF PBZ=150,QUANTITY=''PRESSURE''/']; fprintf(fid,'%s\n',slcf);
fprintf(fid,'%s\n','  '); % blank line

% % find DEVC locations from hill data
% xdevc = [414,651,763,858,920,984,1055,1124,1262];
% ydevc = [1080,1336,1456,1562,1625,1695,1770,1842,1975];
% zdevc = interp2(X,Y,Z_data,xdevc,ydevc)+10;

devc = ['&DEVC ID=''ASW85'', XYZ=414.0,1080.0,18.6, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW50'', XYZ=651.0,1336.0,22.5, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW35'', XYZ=763.0,1456.0,50.0, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW20'', XYZ=858.0,1562.0,85.7, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW10'', XYZ=920.0,1625.0,119.0, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''HT'',    XYZ=984.0,1695.0,134.6, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE10'', XYZ=1055.0,1770.0,117.6, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE20'', XYZ=1124.0,1842.0,94.4, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE40'', XYZ=1262.0,1975.0,51.0, QUANTITY=''VELOCITY''/']; fprintf(fid,'%s\n',devc);

devc = ['&DEVC ID=''ASW85U'', XYZ=414.0,1080.0,18.6, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW50U'', XYZ=651.0,1336.0,22.5, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW35U'', XYZ=763.0,1456.0,50.0, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW20U'', XYZ=858.0,1562.0,85.7, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW10U'', XYZ=920.0,1625.0,119.0, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''HTU'',    XYZ=984.0,1695.0,134.6, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE10U'', XYZ=1055.0,1770.0,117.6, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE20U'', XYZ=1124.0,1842.0,94.4, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE40U'', XYZ=1262.0,1975.0,51.0, QUANTITY=''U-VELOCITY''/']; fprintf(fid,'%s\n',devc);

devc = ['&DEVC ID=''ASW85V'', XYZ=414.0,1080.0,18.6, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW50V'', XYZ=651.0,1336.0,22.5, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW35V'', XYZ=763.0,1456.0,50.0, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW20V'', XYZ=858.0,1562.0,85.7, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ASW10V'', XYZ=920.0,1625.0,119.0, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''HTV'',    XYZ=984.0,1695.0,134.6, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE10V'', XYZ=1055.0,1770.0,117.6, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE20V'', XYZ=1124.0,1842.0,94.4, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);
devc = ['&DEVC ID=''ANE40V'', XYZ=1262.0,1975.0,51.0, QUANTITY=''V-VELOCITY''/']; fprintf(fid,'%s\n',devc);

fprintf(fid,'%s\n','  '); % blank line

% for i=1:N(1)
%     for j=1:N(2)

%         xi = xf(i)+dx/2;
%         yi = yf(j)+dy/2;

%         if xi>=min_x_data & xi<=max_x_data & yi>=min_y_data & yi<=max_y_data

%             zi = interp2(X,Y,Z_data,xi,yi);
%             obst = ['&OBST XB=',num2str(xf(i),'%1.2f'),',',num2str(xf(i)+dx,'%1.2f'),',',...
%                 num2str(yf(j),'%1.2f'),',',num2str(yf(j)+dy,'%1.2f'),',',...
%                 num2str(z_min,'%1.2f'),',',num2str(zi,'%1.2f'),',SURF_ID=''terrain''/'];

%             fprintf(fid,'%s\n',obst);

%         else

%             zi = min_z_data;
%             obst = ['&OBST XB=',num2str(xf(i),'%1.2f'),',',num2str(xf(i)+dx,'%1.2f'),',',...
%                     num2str(yf(j),'%1.2f'),',',num2str(yf(j)+dy,'%1.2f'),',',...
%                     num2str(z_min,'%1.2f'),',',num2str(zi,'%1.2f'),',SURF_ID=''terrain''/'];

%             fprintf(fid,'%s\n',obst);

%         end

%     end
% end

% for GEOM option
ivals = num2str(N(1)+1);
jvals = num2str(N(2)+1);
x1 = num2str(x_min,'%1.2f');
x2 = num2str(x_max,'%1.2f');
y1 = num2str(y_min,'%1.2f');
y2 = num2str(y_max,'%1.2f');

matl = ['&MATL ID=''terrain'', CONDUCTIVITY=1, DENSITY=1000, SPECIFIC_HEAT=1/']; fprintf(fid,'%s\n',matl);

fprintf(fid,'%s\n','  '); % blank line

geom = ['&GEOM ID=''terrain'', MATL_ID=''terrain''']  ;  fprintf(fid,'%s\n',geom);
geom = ['      IJK=',ivals,',',jvals]                 ;  fprintf(fid,'%s\n',geom);
geom = ['      XB=',x1,',',x2,',',y1,',',y2]          ;  fprintf(fid,'%s\n',geom);
geom = ['      ZVALS=']                               ;  fprintf(fid,'%s\n',geom);
geom = ['      ']                                     ;  fprintf(fid,'%s',geom);

n = 0;
for j=1:N(2)+1
    for i=1:N(1)+1

        if xf(i)>min_x_data & xf(i)<max_x_data & yf(j)>min_y_data & yf(j)<max_y_data
            zi = interp2(X,Y,Z_data,xf(i),yf(j));
        else
            zi = min_z_data;
        end

        n = n+1;
        if mod(n,10)==0
            zvals = [num2str(zi,'%1.2f')]; fprintf(fid,'%s\n',zvals);
        else
            zvals = [num2str(zi,'%1.2f')]; fprintf(fid,'%s,',zvals);
        end

    end
end
geom = ['      /']; fprintf(fid,'%s',geom);
n

fprintf(fid,'%s\n','  '); % blank line

tail = ['&TAIL /']; fprintf(fid,'%s\n',tail);

fclose(fid);
