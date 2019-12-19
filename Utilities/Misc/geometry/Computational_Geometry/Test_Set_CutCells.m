close all 
clear all
clc

% Declarations, call set constants:
global MAX_DIM IAXIS JAXIS KAXIS NOD1 NOD2 NOD3 NOD4
global NMESHES MESHES BODINT_PLANE
global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_BDNUM_CRS_AUX IBM_IS_CRS2_AUX
global NM N_GEOMETRY GEOM IBM_INBOUNDARY
global IBM_IDCF IBM_IDCE
global XCELL DXCELL XFACE DXFACE 
global YCELL DYCELL YFACE DYFACE
global ZCELL DZCELL ZFACE DZFACE

global LOW_IND HIGH_IND

global basedir CELLRT

plot_cutedges=false;

[ierr]=SET_CONSTANTS();

% Case name and directory:
addpath ../
basedir = '/Users/mnv/Documents/FIREMODELS_FORK/fds/GEOM_Intersection/';
casename= 'two_spheres';

% Load and plot Geometries:
[GEOM,N_GEOMETRY]=load_geometries(basedir,casename);
%newfig =1;
%[ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig);

% Load and plot Meshes:
[MESHES,NMESHES]=load_meshes(basedir,casename);
newfig=0;
%[ierr]=plot_meshes(NMESHES,MESHES,newfig);

tic
[ierr]=SET_CUTCELLS_3D(basedir,casename,plot_cutedges);
toc


figure
axis equal; box on;
xlabel('X')
ylabel('Y')
zlabel('Z')
NM=1;
% Plot Boundary cut-faces:
for ICF=1:MESHES(NM).N_CUTFACE_MESH
   if(MESHES(NM).CUT_FACE(ICF).STATUS ~= IBM_INBOUNDARY); continue; end
   NFACE=MESHES(NM).CUT_FACE(ICF).NFACE;
   XYZVERT = MESHES(NM).CUT_FACE(ICF).XYZVERT;
   %if(min(XYZVERT(JAXIS,1:MESHES(NM).CUT_FACE(ICF).NVERT)) < 0.); continue; end
   for JCF=1:NFACE
       NELEM  = MESHES(NM).CUT_FACE(ICF).CFELEM(1,JCF);
       CFELEM = MESHES(NM).CUT_FACE(ICF).CFELEM(2:NELEM+1,JCF);
       
       [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'b');
       set(hp,'FaceAlpha',0.3)
   end
   
%    I = MESHES(NM).CUT_FACE(ICF).IJK(IAXIS);
%    J = MESHES(NM).CUT_FACE(ICF).IJK(JAXIS);
%    K = MESHES(NM).CUT_FACE(ICF).IJK(KAXIS);
%    P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
%    P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
%    P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
%    P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
%    P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
%    P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
%    P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
%    P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];
%    
%    FC(1,:) = [ 1 4 8 5 ];
%    FC(2,:) = [ 2 3 7 6 ];
%    FC(3,:) = [ 1 5 6 2 ];
%    FC(4,:) = [ 4 8 7 3 ];
%    FC(5,:) = [ 1 2 3 4 ];
%    FC(6,:) = [ 5 6 7 8 ];
%    
%    for IFC=1:6
%        [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
%        set(hp,'FaceAlpha',0.3)
%    end
end

% a=[];
% for ICF=1:MESHES(NM).N_CUTFACE_MESH
%    if(MESHES(NM).CUT_FACE(ICF).STATUS ~= IBM_INBOUNDARY); continue; end
%    NFACE=MESHES(NM).CUT_FACE(ICF).NFACE;
%    NBD  =length(unique(MESHES(NM).CUT_FACE(ICF).BODTRI(1,1:NFACE)));
%    if( NBD > 1)
%        a =[ a ICF];
%        disp(['ICF=' num2str(ICF)])
%    end
% end

% [ierr]=plot_meshes(NMESHES,MESHES,newfig);


% Plot special cells:
disp(['N_SPCELL=' num2str(MESHES(NM).N_SPCELL)])
for ICELL=1:MESHES(NM).N_SPCELL
    I=MESHES(NM).SPCELL_LIST(IAXIS,ICELL);
    J=MESHES(NM).SPCELL_LIST(JAXIS,ICELL);
    K=MESHES(NM).SPCELL_LIST(KAXIS,ICELL);
    
    P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
    P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
    P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
    P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
    P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
    P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
    P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
    P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];

    FC(1,:) = [ 1 4 8 5 ];
    FC(2,:) = [ 2 3 7 6 ];
    FC(3,:) = [ 1 5 6 2 ];
    FC(4,:) = [ 4 8 7 3 ];
    FC(5,:) = [ 1 2 3 4 ];
    FC(6,:) = [ 5 6 7 8 ];

    for IFC=1:6
    [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
       set(hp,'FaceAlpha',0.3)
    end
end




return