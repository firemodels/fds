function [ierr]=SET_CONSTANTS()

global MAX_DIM IAXIS JAXIS KAXIS NOD1 NOD2 NOD3 NOD4 NGUARD NGUARD_FDS CCGUARD
global GEOMEPS IBM_SOLID IBM_CUTCFE IBM_GASPHASE IBM_UNDEFINED IBM_INBOUNDARY
global IBM_INBOUNDCC IBM_INBOUNDCF
global LOW_IND HIGH_IND
global NODS_WSEL EDG1 EDG2 EDG3

global IBM_N_CRS IBM_MAXCROSS_X2

global IBM_GG IBM_SS IBM_GS IBM_SG IBM_SGG

global IBM_VGSC IBM_NVVARS
global IBM_EGSC IBM_IDCE IBM_ECRS IBM_NEVARS
global IBM_FGSC IBM_IDCF IBM_FFNF IBM_IDRA IBM_NFVARS
global IBM_CGSC IBM_IDCC IBM_UNKZ IBM_UNKH IBM_NCVARS
global IBM_FTYPE_RGGAS IBM_FTYPE_CFGAS IBM_FTYPE_CFINB IBM_FTYPE_SVERT ...
       IBM_FTYPE_RCGAS IBM_FTYPE_CCGAS
global FCELL
global IBM_MAX_WSTRIANG_SEG DELTA_SEGBIN DELTA_VERT
global IBM_DELTA_NBCROSS IBM_MAXVERTS_FACE IBM_MAXCEELEM_FACE

global IBM_MAXVERTS_CELL IBM_NPARAM_CCFACE

global DELTA_EDGE DELTA_FACE DELTA_CELL

ierr=1;

MAX_DIM=3;
IAXIS  =1;
JAXIS  =2;
KAXIS  =3;
NOD1   =1;
NOD2   =2;
NOD3   =3;
NOD4   =4;
EDG1   =1;
EDG2   =2;
EDG3   =3;
FCELL  =1;
NODS_WSEL=3;
IBM_MAX_WSTRIANG_SEG=2;
DELTA_SEGBIN = 50;
IBM_DELTA_NBCROSS=20;
DELTA_VERT   = 24;
DELTA_EDGE   = 24;
DELTA_FACE   = 24; 
DELTA_CELL   =  5;

LOW_IND = 1;
HIGH_IND= 2;

NGUARD_FDS   = 1;
NGUARD       = 6;
CCGUARD      = NGUARD-2;

GEOMEPS = 1.e-12;
IBM_INBOUNDARY=2;
IBM_SOLID    = 1;
IBM_CUTCFE   = 0;
IBM_GASPHASE =-1;
IBM_INBOUNDCF=-2;
IBM_INBOUNDCC=-3;
IBM_UNDEFINED=-11;

IBM_GG =  1; % Gas - Gas intersection.
IBM_SS =  3; % Solid - Solid intersection.
IBM_GS = -1; % Gas to Solid intersection (as coordinate xi increases).
IBM_SG =  5; % Solid to Gas intersection (as coordinate xi increases).
IBM_SGG= IBM_GG; % Single point GG intersection. Might not be needed.

% Constants used to identify variables on Eulerian grid arrays:
% Vertex centered variables:
IBM_VGSC   = 1; % Type of vertex media, IBM_GASPHASE or IBM_SOLID.
IBM_NVVARS = 1; % Number of vertex variables in MESHES(N)%IBM_VERTVAR.

% Cartesian Edge centered variables:
IBM_EGSC   = 1; % Edge media type: IBM_GASPHASE, IBM_SOLID or IBM_CUTCFE.
IBM_IDCE   = 2; % MESHES(N)%CUT_EDGE data struct entry index location.
IBM_ECRS   = 3; % MESHES(N)%EDGE_CROSS data struct entry index location.
IBM_NEVARS = 3; % Number of edge variables in MESHES(N)%ECVAR.

% Cartesian Face centered variables:
IBM_FGSC   = 1; % Face media type: IBM_GASPHASE, IBM_SOLID or IBM_CUTCFE.
%IBM_IDCE   = 2 % MESHES(N)%CUT_EDGE data struct entry index location,
                                      % IBM_INBOUNDCF type.
IBM_IDCF   = 3; % MESHES(N)%CUT_FACE data struct entry index location,
                                      % IBM_INBOUNDCF type cut-faces.
IBM_FFNF   = 4; % Flag that defines if face is to be IB forced or not.
IBM_IDRA   = 5; % Integer ICF that defines RAD_FACE(ICF) position for radiation FVM.
IBM_NFVARS = 5; % Number of face variables in MESHES(N)%FCVAR.

% Cartesian Cell centered variables:
IBM_CGSC   = 1; % Face media type: IBM_GASPHASE, IBM_SOLID or IBM_CUTCFE.
%IBM_IDCE   = 2 % MESHES(N)%CUT_EDGE data struct entry index location,
                                      % cut edges in Cartesian cell.
%IBM_IDCF   = 3 % MESHES(N)%CUT_FACE data struct entry index location,
                                      % IBM_INBOUNDARY type cut-faces (and CFACEs) in Cartesian cell.
IBM_IDCC   = 4; % MESHES(N)%CUT_CELL data struct entry index location,
                                      % cut-cells in Cartesian cell.
IBM_UNKZ   = 5; % Scalar indexing.
IBM_UNKH   = 6; % Pressure indexing.
IBM_NCVARS = 6; % Number of cell variables in MESHES(N)%CCVAR.

% Cut-faces types in FACE_LIST of CUT_CELL:
IBM_FTYPE_RGGAS = 0; % This face of a cut-cell is a regular GASPHASE face.
IBM_FTYPE_CFGAS = 1; % A GASPHASE cut-face or cell.
IBM_FTYPE_CFINB = 2; % An INBOUNDARY cut-face.
IBM_FTYPE_SVERT = 3; % A SOLID Vertex.
% Extra tagging parameters, for RESTRICT_EP:
IBM_FTYPE_RCGAS = 4; % A Face between a regular cell and a cut-cell.
IBM_FTYPE_CCGAS = 5; % A regular gas cut-cell.

IBM_MAXCROSS_X2 = 3072;
IBM_MAXVERTS_FACE=3072;
IBM_MAXCEELEM_FACE=3072;
IBM_N_CRS = 0;

IBM_MAXVERTS_CELL   =3072;
IBM_NPARAM_CCFACE   =   5; 

ierr=0;

return