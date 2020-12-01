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

NM=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM=NM+1;
        M(NM).xs = XB(1) + i*DX;
        M(NM).xf = XB(2) + i*DX;
        M(NM).ys = XB(3) + j*DY;
        M(NM).yf = XB(4) + j*DY;
        M(NM).x = M(NM).xs:dx:M(NM).xf;
        M(NM).y = M(NM).ys:dy:M(NM).yf;
        [M1(NM).X,M1(NM).Y] = meshgrid(M(NM).x,M(NM).y);
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

NM=0;
for j=J_LOWER:J_UPPER
    for i=I_LOWER:I_UPPER
        NM=NM+1;
        M(NM).xs = XB(1) + i*DX;
        M(NM).xf = XB(2) + i*DX;
        M(NM).ys = XB(3) + j*DY;
        M(NM).yf = XB(4) + j*DY;
        M(NM).x = M(NM).xs:dx:M(NM).xf;
        M(NM).y = M(NM).ys:dy:M(NM).yf;
        [M2(NM).X,M2(NM).Y] = meshgrid(M(NM).x,M(NM).y);
    end
end