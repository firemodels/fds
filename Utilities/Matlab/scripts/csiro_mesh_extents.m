% McDermott
% 12-1-2020
% csiro_mesh_extents.m

% from Case_F19 input file
% &MESH IJK=80,80,40, XB=-20.0,20.0,0.0,40.0,0.0,20.0, MULT_ID='mesh' /
% &MULT ID='mesh', DX=40., DY=40., DZ=20., I_LOWER=0, I_UPPER=5, J_LOWER=-3, J_UPPER=2, K_UPPER=0 /

XB=[-20.0,20.0,0.0,40.0,0.0,20.0];
NX=80;
NY=80;
DX=40;
DY=40;
dx=DX/NX;
dy=DY/NY;
I_LOWER=0;
I_UPPER=5;
J_LOWER=-3;
J_UPPER=2;

NM1=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM1=NM1+1;
        M(NM1).xs = XB(1) + i*DX;
        M(NM1).xf = XB(2) + i*DX;
        M(NM1).ys = XB(3) + j*DY;
        M(NM1).yf = XB(4) + j*DY;
        M(NM1).x = M(NM1).xs:dx:M(NM1).xf;
        M(NM1).y = M(NM1).ys:dy:M(NM1).yf;
        [M1(NM1).X,M1(NM1).Y] = meshgrid(M(NM1).x,M(NM1).y);
    end
end

% from Case_C064 input file
% &MESH IJK=40,40,40, XB=-10.0,10.0,0.0,20.0,0.0,20.0, MULT_ID='mesh' /
% &MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=5, J_LOWER=-3, J_UPPER=2, K_UPPER=0 /

XB=[-10.0,10.0,0.0,20.0,0.0,20.0];
NX=40;
NY=40;
DX=20;
DY=20;
dx=DX/NX;
dy=DY/NY;
I_LOWER=0;
I_UPPER=5;
J_LOWER=-3;
J_UPPER=2;

NM2=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM2=NM2+1;
        M(NM2).xs = XB(1) + i*DX;
        M(NM2).xf = XB(2) + i*DX;
        M(NM2).ys = XB(3) + j*DY;
        M(NM2).yf = XB(4) + j*DY;
        M(NM2).x = M(NM2).xs:dx:M(NM2).xf;
        M(NM2).y = M(NM2).ys:dy:M(NM2).yf;
        [M2(NM2).X,M2(NM2).Y] = meshgrid(M(NM2).x,M(NM2).y);
    end
end

% from Case_F19_fine input file
% &MESH IJK=80,80,80, XB=-20.0,0.0,0.0,20.0,0.0,20.0, MULT_ID='mesh' /
% &MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=11, J_LOWER=-6, J_UPPER=5, K_UPPER=0 /

XB=[-20.0,0.0,0.0,20.0,0.0,20.0];
NX=80;
NY=80;
DX=20;
DY=20;
dx=DX/NX;
dy=DY/NY;
I_LOWER=0;
I_UPPER=11;
J_LOWER=-6;
J_UPPER=5;

NM3=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM3=NM3+1;
        M(NM3).xs = XB(1) + i*DX;
        M(NM3).xf = XB(2) + i*DX;
        M(NM3).ys = XB(3) + j*DY;
        M(NM3).yf = XB(4) + j*DY;
        M(NM3).x = M(NM3).xs:dx:M(NM3).xf;
        M(NM3).y = M(NM3).ys:dy:M(NM3).yf;
        [M3(NM3).X,M3(NM3).Y] = meshgrid(M(NM3).x,M(NM3).y);
    end
end

% from Case_C064_fine input file
% &MESH IJK=80,80,80, XB=-10.0,10.0,0.0,20.0,0.0,20.0, MULT_ID='mesh' /
% &MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=5, J_LOWER=-3, J_UPPER=2, K_UPPER=0 /

XB=[-10.0,10.0,0.0,20.0,0.0,20.0];
NX=80;
NY=80;
DX=20;
DY=20;
dx=DX/NX;
dy=DY/NY;
I_LOWER=0;
I_UPPER=5;
J_LOWER=-3;
J_UPPER=2;

NM4=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM4=NM4+1;
        M(NM4).xs = XB(1) + i*DX;
        M(NM4).xf = XB(2) + i*DX;
        M(NM4).ys = XB(3) + j*DY;
        M(NM4).yf = XB(4) + j*DY;
        M(NM4).x = M(NM4).xs:dx:M(NM4).xf;
        M(NM4).y = M(NM4).ys:dy:M(NM4).yf;
        [M4(NM4).X,M4(NM4).Y] = meshgrid(M(NM4).x,M(NM4).y);
    end
end
