% Test to load closed terrain and strip lower set of triangles. The result
% is the terrain surface trinagulation only.
%
% -------------------------------------------------------------------------
close all
clear all

IAXIS=1; JAXIS=2; KAXIS=3;
NOD1 =1; NOD2=2;  NOD3=3; 

basedir='/Volumes/mnv/FIREMODELS_FORK/fds/Meyrin_Cases/ISOLDE/Grid_5m/';
casename='02_Isolde_234msh_5m_5ms';
geomname='MNS';

[READFLG,GEOM]=load_bingeom(basedir,casename,geomname);

% First figure:
figure
trisurf(GEOM.WSELEM,GEOM.XYZ(:,IAXIS),GEOM.XYZ(:,JAXIS),GEOM.XYZ(:,KAXIS),'Edgecolor','interp')
axis equal


% Strip Triangles:
% First define nodes INOD which lay below Z_THRES in NODE_ID(INOD)=0: 
Z_THRES = 15;
NODE_ID  = zeros(1,GEOM.N_VERTS);
NODE_ID2 = zeros(1,GEOM.N_VERTS);
[I]=find(GEOM.XYZ(:,KAXIS) > Z_THRES);
NODE_ID(I)=1;

% Renumber nodes:
COUNT=0;
for INOD=1:GEOM.N_VERTS
   if(NODE_ID(INOD) == 0); continue; end
   COUNT=COUNT+1;
   NODE_ID(INOD)  = COUNT;
   NODE_ID2(COUNT)= INOD;
end

% Then Add nodes to output geometry GEOM2: 
GEOM2.N_VERTS=COUNT;
GEOM2.XYZ = GEOM.XYZ(NODE_ID2(1:COUNT),IAXIS:KAXIS);

% Now strip Renumber output faces:
WSELEM= zeros(GEOM.N_FACES,NOD3);
SURFS = zeros(1,GEOM.N_FACES);
COUNT=0;
for IFC=1:GEOM.N_FACES
   IWSELEM = GEOM.WSELEM(IFC,NOD1:NOD3);
   NODF_ID= NODE_ID(IWSELEM);
   if(any(NODF_ID==0)); continue; end
   COUNT = COUNT + 1;
   WSELEM(COUNT,NOD1:NOD3) = NODF_ID;
   SURFS(COUNT) = GEOM.SURFS(IFC);
end
GEOM2.N_FACES=COUNT;
GEOM2.WSELEM=WSELEM(1:COUNT,NOD1:NOD3);
GEOM2.SURFS=SURFS(1:COUNT);
GEOM2.N_SURF_ID=GEOM.N_SURF_ID;
GEOM2.N_VOLUS=0;

GEOM2.VERTS = zeros(1,3*GEOM2.N_VERTS);
for INOD=1:GEOM2.N_VERTS
    GEOM2.VERTS(KAXIS*(INOD-1)+1:KAXIS*INOD) = GEOM2.XYZ(INOD,IAXIS:KAXIS);
end
GEOM2.FACES = zeros(1,3*GEOM2.N_FACES);
for IFC=1:GEOM2.N_FACES
    GEOM2.FACES(NOD3*(IFC-1)+1:NOD3*IFC) = GEOM2.WSELEM(IFC,NOD1:NOD3);
end

% Plot stripped terrain:
figure
trisurf(GEOM2.WSELEM,GEOM2.XYZ(:,IAXIS),GEOM2.XYZ(:,JAXIS),GEOM2.XYZ(:,KAXIS),'Edgecolor','interp')
axis equal


% Finally write casename_geomname_OUT.bingeom:
[WRITEFLG]=write_bingeom(GEOM2,basedir,casename,geomname);

return