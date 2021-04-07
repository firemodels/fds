% McDermott
% 13 Feb 2021
% pan_with_lip.m
% conversion of pan_with_lip.f90 (mcgratta) to matlab
%
% Estimate area of intersection of circle with origin (X0,Y0), radius, RAD, and rectangle (X1,Y1,X2,Y2)
% All dimensions in meters

X0 = 0.;
Y0 = 0.;
Z1 = -0.2;             % NIST 1 m pan depth (15 cm), use 0.2 m to extend out of domain for HT3D
Z2 = 0.;
RAD = 0.503;           % NIST 1 m outside diameter
PAN_THICKNESS = 0.01;  % NIST 1 m thickness (0.17 cm), use max(THICKNESS,DX)
LIP_HEIGHT = 0.01;     % NIST 1 m height (1.0 cm), use max(HEIGHT,DZ)

IBAR = 120; % dx=4 cm -->IBAR=30; dx=2 cm --> IBAR=60; dx=1 cm --> IBAR=120
JBAR = 120;
XS = -0.6;
XF = 0.6;
YS = -0.6;
YF = 0.6;

RAD2 = RAD + PAN_THICKNESS;
DX = (XF-XS)/double(IBAR);
DY = (YF-YS)/double(JBAR);
AREA = 0.;

fid = fopen('nist_1m_pan_dx1cm.txt','w');

for J=1:JBAR
    Y = YS + DY*(J-0.5);
    for I=1:IBAR
        X = XS + DX*(I-0.5);
        if ((X-X0)^2+(Y-Y0)^2)<RAD^2
            fprintf(fid,'%s %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %s \n', ...
                "&OBST XB=",X-0.5*DX,X+0.5*DX,Y-0.5*DY,Y+0.5*DY,Z1,Z2, " HT3D=T, SURF_ID='POOL' /");
            AREA = AREA + DX*DY;
        elseif ((X-X0)^2+(Y-Y0)^2)<RAD2^2
            fprintf(fid,'%s %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %s \n', ...
                "&OBST XB=",X-0.5*DX,X+0.5*DX,Y-0.5*DY,Y+0.5*DY,Z1,Z2+LIP_HEIGHT, " HT3D=T, SURF_ID='PAN' /");
        end
    end
end

disp(['Area of pan is ',num2str(AREA),' m2'])

fclose(fid);
