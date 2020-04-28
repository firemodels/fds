function [ierr]=SET_CUTCELLS_3D(basedir,casename,plot_cutedges)

global NMESHES MESHES NGUARD CCGUARD GEOMEPS IAXIS JAXIS KAXIS NOD1 NOD2 MAX_DIM
global X1FACE X2FACE X3FACE DX1FACE DX2FACE DX3FACE DX2CELL DX3CELL
global X2LO X2HI X3LO X3HI
global X2CELL X3CELL X2LO_CELL X2HI_CELL X3LO_CELL X3HI_CELL
global ILO_FACE IHI_FACE ILO_CELL IHI_CELL
global JLO_FACE JHI_FACE JLO_CELL JHI_CELL
global KLO_FACE KHI_FACE KLO_CELL KHI_CELL
global XCELL DXCELL XFACE DXFACE
global YCELL DYCELL YFACE DYFACE
global ZCELL DZCELL ZFACE DZFACE
global BODINT_PLANE
global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS IBM_INBOUNDARY IBM_GASPHASE IBM_CGSC
global X1NOC X2NOC X3NOC TRANS
global X1LO_CELL X1HI_CELL
global FACERT CELLRT
global N_GEOMETRY GEOM

global XIAXIS XJAXIS XKAXIS

ierr=1;

for NM=1:NMESHES

   disp(['SET_CUTCELLS_3D : Processing mesh ' num2str(NM)])

   IBAR = MESHES(NM).IBAR; JBAR = MESHES(NM).JBAR; KBAR = MESHES(NM).KBAR;

   % Mesh sizes:
   NXB=IBAR;   NYB=JBAR;   NZB=KBAR;

   % X direction bounds:
   ILO_FACE = MESHES(NM).ILO_FACE;  % Low mesh boundary face index.
   IHI_FACE = MESHES(NM).IHI_FACE;  % High mesh boundary face index.
   ILO_CELL = MESHES(NM).ILO_CELL;  % First internal cell index. See notes.
   IHI_CELL = MESHES(NM).IHI_CELL;  % Last internal cell index.
   ISTR     = MESHES(NM).ISTR;      % Allocation start x arrays.
   IEND     = MESHES(NM).IEND;      % Allocation end x arrays.

   % Y direction bounds:
   JLO_FACE = MESHES(NM).JLO_FACE;  % Low mesh boundary face index.
   JHI_FACE = MESHES(NM).JHI_FACE;  % High mesh boundary face index.
   JLO_CELL = MESHES(NM).JLO_CELL;  % First internal cell index. See notes.
   JHI_CELL = MESHES(NM).JHI_CELL;  % Last internal cell index.
   JSTR     = MESHES(NM).JSTR;      % Allocation start y arrays.
   JEND     = MESHES(NM).JEND;      % Allocation end y arrays.

   % Z direction bounds:
   KLO_FACE = MESHES(NM).KLO_FACE;  % Low mesh boundary face index.
   KHI_FACE = MESHES(NM).KHI_FACE;  % High mesh boundary face index.
   KLO_CELL = MESHES(NM).KLO_CELL;  % First internal cell index. See notes.
   KHI_CELL = MESHES(NM).KHI_CELL;  % Last internal cell index.
   KSTR     = MESHES(NM).KSTR;      % Allocation start z arrays.
   KEND     = MESHES(NM).KEND;      % Allocation end z arrays.

   % Define grid arrays for this mesh:
   % Populate position and cell size arrays: Uniform grid implementation.
   % X direction:
   DXCELL=zeros(1,IEND); DXCELL(ILO_CELL-1:IHI_CELL+1)=MESHES(NM).DX(ILO_CELL-1:IHI_CELL+1);
   for IGC=2:NGUARD
      DXCELL(ILO_CELL-IGC)=DXCELL(ILO_CELL-IGC+1);
      DXCELL(IHI_CELL+IGC)=DXCELL(IHI_CELL+IGC-1);
   end
   DXFACE=zeros(1,IEND); DXFACE(ILO_FACE:IHI_FACE)=MESHES(NM).DXN(ILO_FACE:IHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      DXFACE(ILO_FACE-IGC)=DXFACE(ILO_FACE-IGC+1);
      DXFACE(IHI_FACE+IGC)=DXFACE(ILO_FACE+IGC-1);
   end
   XCELL=zeros(1,IEND);  XCELL = 1./GEOMEPS; % Initialize huge.
   XCELL(ILO_CELL-1:IHI_CELL+1)=MESHES(NM).XC(ILO_CELL-1:IHI_CELL+1);
   for IGC=2:NGUARD
      XCELL(ILO_CELL-IGC)=XCELL(ILO_CELL-IGC+1)-DXFACE(ILO_FACE-IGC+1);
      XCELL(IHI_CELL+IGC)=XCELL(IHI_CELL+IGC-1)+DXFACE(IHI_FACE+IGC-1);
   end
   XFACE=zeros(1,IEND);  XFACE = 1./GEOMEPS; % Initialize huge.
   XFACE(ILO_FACE:IHI_FACE)=MESHES(NM).X(ILO_FACE:IHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      XFACE(ILO_FACE-IGC)=XFACE(ILO_FACE-IGC+1)-DXCELL(ILO_CELL-IGC);
      XFACE(IHI_FACE+IGC)=XFACE(IHI_FACE+IGC-1)+DXCELL(IHI_CELL+IGC);
   end

   % Y direction:
   DYCELL=zeros(1,JEND); DYCELL(JLO_CELL-1:JHI_CELL+1)= MESHES(NM).DY(JLO_CELL-1:JHI_CELL+1);
   for IGC=2:NGUARD
      DYCELL(JLO_CELL-IGC)=DYCELL(JLO_CELL-IGC+1);
      DYCELL(JHI_CELL+IGC)=DYCELL(JHI_CELL+IGC-1);
   end
   DYFACE=zeros(1,JEND); DYFACE(JLO_FACE:JHI_FACE)= MESHES(NM).DYN(JLO_FACE:JHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      DYFACE(JLO_FACE-IGC)=DYFACE(JLO_FACE-IGC+1);
      DYFACE(JHI_FACE+IGC)=DYFACE(JHI_FACE+IGC-1);
   end
   YCELL=zeros(1,JEND);  YCELL = 1./GEOMEPS; % Initialize huge.
   YCELL(JLO_CELL-1:JHI_CELL+1) = MESHES(NM).YC(JLO_CELL-1:JHI_CELL+1);
   for IGC=2:NGUARD
      YCELL(JLO_CELL-IGC)=YCELL(JLO_CELL-IGC+1)-DYFACE(JLO_FACE-IGC+1);
      YCELL(JHI_CELL+IGC)=YCELL(JHI_CELL+IGC-1)+DYFACE(JHI_FACE+IGC-1);
   end
   YFACE=zeros(1,JEND);  YFACE = 1./GEOMEPS; % Initialize huge.
   YFACE(JLO_FACE:JHI_FACE) = MESHES(NM).Y(JLO_FACE:JHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      YFACE(JLO_FACE-IGC)=YFACE(JLO_FACE-IGC+1)-DYCELL(JLO_CELL-IGC);
      YFACE(JHI_FACE+IGC)=YFACE(JHI_FACE+IGC-1)+DYCELL(JHI_CELL+IGC);
   end

   % Z direction:
   DZCELL=zeros(1,KEND); DZCELL(KLO_CELL-1:KHI_CELL+1)=MESHES(NM).DZ(KLO_CELL-1:KHI_CELL+1);
   for IGC=2:NGUARD
      DZCELL(KLO_CELL-IGC)=DZCELL(KLO_CELL-IGC+1);
      DZCELL(KHI_CELL+IGC)=DZCELL(KHI_CELL+IGC-1);
   end
   DZFACE=zeros(1,KEND); DZFACE(KLO_FACE:KHI_FACE)=MESHES(NM).DZN(KLO_FACE:KHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      DZFACE(KLO_FACE-IGC)=DZFACE(KLO_FACE-IGC+1);
      DZFACE(KHI_FACE+IGC)=DZFACE(KHI_FACE+IGC-1);
   end
   ZCELL=zeros(1,KEND);  ZCELL = 1./GEOMEPS; % Initialize huge.
   ZCELL(KLO_CELL-1:KHI_CELL+1)=MESHES(NM).ZC(KLO_CELL-1:KHI_CELL+1);
   for IGC=2:NGUARD
      ZCELL(KLO_CELL-IGC)=ZCELL(KLO_CELL-IGC+1)-DZFACE(KLO_FACE-IGC+1);
      ZCELL(KHI_CELL+IGC)=ZCELL(KHI_CELL+IGC-1)+DZFACE(KHI_FACE+IGC-1);
   end
   ZFACE=zeros(1,KEND);  ZFACE = 1./GEOMEPS; % Initialize huge.
   ZFACE(KLO_FACE:KHI_FACE) = MESHES(NM).Z(KLO_FACE:KHI_FACE);
   for IGC=1:NGUARD-1 % OJO
      ZFACE(KLO_FACE-IGC)=ZFACE(KLO_FACE-IGC+1)-DZCELL(KLO_CELL-IGC);
      ZFACE(KHI_FACE+IGC)=ZFACE(KHI_FACE+IGC-1)+DZCELL(KHI_CELL+IGC);
   end

   CELLRT=zeros(IEND,JEND,KEND); % false

   % Here Allocation of crossings, cut-edges, faces and cells..
   disp(['Calculation of intersections and cut-edges by grid plane ...'])
   tic
   for X1AXIS=IAXIS:KAXIS

       switch(X1AXIS)
       case(IAXIS)

          PLNORMAL = [ 1., 0., 0.];
          ILO = ILO_FACE-CCGUARD;  IHI = IHI_FACE+CCGUARD;
          JLO = JLO_FACE;  JHI = JLO_FACE;
          KLO = KLO_FACE;  KHI = KLO_FACE;

          % x2, x3 axes parameters:
          X2AXIS = JAXIS; X2LO = JLO_FACE-CCGUARD; X2HI = JHI_FACE+CCGUARD;
          X3AXIS = KAXIS; X3LO = KLO_FACE-CCGUARD; X3HI = KHI_FACE+CCGUARD;

          % location in I,J,K of x2,x2,x3 axes:
          XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS;

          % Face coordinates in x1,x2,x3 axes:
          X1FACE = zeros(1,IEND); DX1FACE = zeros(1,IEND);
          X1FACE = XFACE; DX1FACE = DXFACE;
          X2FACE = zeros(1,JEND); DX2FACE = zeros(1,JEND);
          X2FACE = YFACE; DX2FACE = DYFACE;
          X3FACE = zeros(1,KEND); DX3FACE = zeros(1,KEND);
          X3FACE = ZFACE; DX3FACE = DZFACE;

          X1LO_CELL = ILO_CELL-CCGUARD; X1HI_CELL = IHI_CELL+CCGUARD;

          % x2 cell center parameters:
          X2LO_CELL = JLO_CELL-CCGUARD; X2HI_CELL = JHI_CELL+CCGUARD;
          X2CELL = zeros(1,JEND); DX2CELL=zeros(1,JEND);
          X2CELL = YCELL; DX2CELL = DYCELL;

          % x3 cell center parameters:
          X3LO_CELL = KLO_CELL-CCGUARD; X3HI_CELL = KHI_CELL+CCGUARD;
          X3CELL = zeros(1,KEND); DX3CELL = zeros(1,KEND);
          X3CELL = ZCELL; DX3CELL = DZCELL;

       case(JAXIS)

          PLNORMAL = [ 0., 1., 0.];
          ILO = ILO_FACE;  IHI = ILO_FACE;
          JLO = JLO_FACE-CCGUARD;  JHI = JHI_FACE+CCGUARD;
          KLO = KLO_FACE;  KHI = KLO_FACE;

          % x2, x3 axes parameters:
          X2AXIS = KAXIS; X2LO = KLO_FACE-CCGUARD; X2HI = KHI_FACE+CCGUARD;
          X3AXIS = IAXIS; X3LO = ILO_FACE-CCGUARD; X3HI = IHI_FACE+CCGUARD;

          % location in I,J,K of x2,x2,x3 axes:
          XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS;

          % Face coordinates in x1,x2,x3 axes:
          X1FACE = zeros(1,JEND); DX1FACE = zeros(1,JEND);
          X1FACE = YFACE; DX1FACE = DYFACE;
          X2FACE = zeros(1,KEND); DX2FACE = zeros(1,KEND);
          X2FACE = ZFACE; DX2FACE = DZFACE;
          X3FACE = zeros(1,IEND); DX3FACE = zeros(1,IEND);
          X3FACE = XFACE; DX3FACE = DXFACE;

          X1LO_CELL = JLO_CELL-CCGUARD; X1HI_CELL = JHI_CELL+CCGUARD;

          % x2 cell center parameters:
          X2LO_CELL = KLO_CELL-CCGUARD; X2HI_CELL = KHI_CELL+CCGUARD;
          X2CELL = zeros(1,KEND); DX2CELL = zeros(1,KEND);
          X2CELL = ZCELL; DX2CELL = DZCELL;

          % x3 cell center parameters:
          X3LO_CELL = ILO_CELL-CCGUARD; X3HI_CELL = IHI_CELL+CCGUARD;
          X3CELL = zeros(1,IEND); DX3CELL = zeros(1,IEND);
          X3CELL = XCELL; DX3CELL = DXCELL;

       case(KAXIS)

          PLNORMAL = [ 0., 0., 1.];
          ILO = ILO_FACE;  IHI = ILO_FACE;
          JLO = JLO_FACE;  JHI = JLO_FACE;
          KLO = KLO_FACE-CCGUARD;  KHI = KHI_FACE+CCGUARD;

          % x2, x3 axes parameters:
          X2AXIS = IAXIS; X2LO = ILO_FACE-CCGUARD; X2HI = IHI_FACE+CCGUARD;
          X3AXIS = JAXIS; X3LO = JLO_FACE-CCGUARD; X3HI = JHI_FACE+CCGUARD;

          % location in I,J,K of x2,x2,x3 axes:
          XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS;

          % Face coordinates in x1,x2,x3 axes:
          X1FACE = zeros(1,KEND); DX1FACE = zeros(1,KEND);
          X1FACE = ZFACE; DX1FACE = DZFACE;
          X2FACE = zeros(1,IEND); DX2FACE = zeros(1,IEND);
          X2FACE = XFACE; DX2FACE = DXFACE;
          X3FACE = zeros(1,JEND); DX3FACE = zeros(1,JEND);
          X3FACE = YFACE; DX3FACE = DYFACE;

          X1LO_CELL = KLO_CELL-CCGUARD; X1HI_CELL = KHI_CELL+CCGUARD;

          % x2 cell center parameters:
          X2LO_CELL = ILO_CELL-CCGUARD; X2HI_CELL = IHI_CELL+CCGUARD;
          X2CELL = zeros(1,IEND); DX2CELL = zeros(1,IEND);
          X2CELL = XCELL; DX2CELL = DXCELL;

          % x3 cell center parameters:
          X3LO_CELL = JLO_CELL-CCGUARD; X3HI_CELL = JHI_CELL+CCGUARD;
          X3CELL = zeros(1,JEND); DX3CELL = zeros(1,JEND);
          X3CELL = YCELL; DX3CELL = DYCELL;

       end

       % Variable that states if raytracing is necessary to define segments
       % status in a cartesian face.
       FACERT = zeros(X2HI_CELL,X3HI_CELL);


       % Stretched grid vars:
       X1NOC=TRANS(NM).NOC(X1AXIS);
       X2NOC=TRANS(NM).NOC(X2AXIS);
       X3NOC=TRANS(NM).NOC(X3AXIS);

       % Loop Coordinate Planes:
       for K=KLO:KHI
          for J=JLO:JHI
             for I=ILO:IHI

               % Which Plane?
               INDX1(IAXIS:KAXIS) = [ I, J, K ];
               X1PLN = X1FACE(INDX1(X1AXIS));

               % Get intersection of body on plane defined by X1PLN, normal to X1AXIS:
               DX2_MIN = min(DX2CELL(X2LO_CELL:X2HI_CELL));
               DX3_MIN = min(DX3CELL(X3LO_CELL:X3HI_CELL));
               TRI_ONPLANE_ONLY = false;
               RAYTRACE_X2_ONLY = false;
               FACERT(:,:) = 0;

               [ierr,BODINT_PLANE]=GET_BODINT_PLANE(X1AXIS,X1PLN,INDX1(X1AXIS),PLNORMAL,X2AXIS,...
                                   X3AXIS,DX2_MIN,DX3_MIN,X2LO,X2HI,X3LO,X3HI,X2FACE,X3FACE,TRI_ONPLANE_ONLY,RAYTRACE_X2_ONLY);

               % Test that there is an intersection:
               if ((BODINT_PLANE.NSGLS+BODINT_PLANE.NSEGS+BODINT_PLANE.NTRIS) == 0); continue; end

               % Drop if node locations outside block plane area:
               if ((X2FACE(X2LO)-max(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))) > GEOMEPS); continue; end
               if ((min(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))-X2FACE(X2HI)) > GEOMEPS); continue; end
               if ((X3FACE(X3LO)-max(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))) > GEOMEPS); continue; end
               if ((min(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))-X3FACE(X3HI)) > GEOMEPS); continue; end

               if (plot_cutedges && X1AXIS==JAXIS); CUT_EDGE_START = MESHES(NM).N_CUTEDGE_MESH; end

               % For plane normal to X1AXIS, shoot rays along X2AXIS on all X3AXIS gridline
               % locations, get intersection data: Loop x3 axis locations
               for KK=X3LO:X3HI

                  % x3 location of ray along x2, on the x2-x3 plane:
                  X3RAY = X3FACE(KK);

                  % Intersections along x2 for X3RAY x3 location:
                  [ierr] = GET_X2_INTERSECTIONS(X1AXIS,X2AXIS,X3AXIS,X3RAY,X1PLN);

                  % Drop x2 ray if all intersections are outside of the MESH block domain:
                  if (IBM_N_CRS > 0)
                     if ((X2FACE(X2LO)-IBM_SVAR_CRS(IBM_N_CRS)) > GEOMEPS)
                        continue
                     elseif (IBM_SVAR_CRS(1)-X2FACE(X2HI) > GEOMEPS)
                        continue
                     end
                  end

                  % Now for this ray, set vertex types in MESHES(NM).VERTVAR(:,:,:,IBM_VGSC):
                  [ierr]=GET_X2_VERTVAR(X1AXIS,X2LO,X2HI,NM,I,KK);

                  % Now define Crossings on Cartesian Edges and Body segments:
                  % Cartesian cut-edges:
                  [ierr]=GET_CARTEDGE_CUTEDGES(X1AXIS,X2AXIS,X3AXIS,XIAXIS,XJAXIS,XKAXIS, ...
                                               NM,X2LO_CELL,X2HI_CELL,INDX1,KK);

                  % Set segment crossings:
                  % This data is defined by plane, add to current:
                  % - BODINT_PLANE : Data structure with information for crossings on
                  %                      body segments.
                  %                   % NBCROSS(1:NSEGS)        = Number of crossings
                  %                                               on the segment.
                  %                   % SVAR(1:NBCROSS,1:NSEGS) = distance from node 1
                  %                                               along the segment.
                  [ierr]=GET_BODX2_INTERSECTIONS(X2AXIS,X3AXIS,X3RAY);

               end % KK - x3 gridlines.

               % Now for segments not aligned with x3, define
               % intersections with grid line vertices:
               [ierr]=GET_BODX3_INTERSECTIONS(X2AXIS,X3AXIS,X2LO,X2HI);

               if (plot_cutedges && X1AXIS==JAXIS)
                   [ierr]=plot_bodintplane(BODINT_PLANE,false);
                   grid off
                   box on
                   title(['X1AXIS=' num2str(X1AXIS) ', X1PLN=' num2str(X1PLN) ', SLICE=' num2str(INDX1(X1AXIS))])
                   xlabel('X2')
                   ylabel('X3')
                   axis([X2FACE(X2LO) X2FACE(X2HI) X3FACE(X3LO) X3FACE(X3HI)])
                   %CUT_EDGE_START = MESHES(NM).N_CUTEDGE_MESH;
               end

               % After these loops all segments should contain points from Node1,
               % cross 1, cross 2, ..., Node2, in ascending sbod order.
               % Time to generate the body IBM_INBOUNDARY edges on faces and add
               % to MESHES(NM).CUT_EDGE:
               [ierr]=GET_CARTFACE_CUTEDGES(X1AXIS,X2AXIS,X3AXIS,            ...
                                    XIAXIS,XJAXIS,XKAXIS,NM,                 ...
                                    X2LO,X2HI,X3LO,X3HI,X2LO_CELL,X2HI_CELL, ...
                                    X3LO_CELL,X3HI_CELL,INDX1,X1PLN);

               if (plot_cutedges && X1AXIS==JAXIS)
                  CUT_EDGE_END = MESHES(NM).N_CUTEDGE_MESH;
                  for KK=X3LO:X3HI
                      plot([X2FACE(X2LO) X2FACE(X2HI)],[X3FACE(KK) X3FACE(KK)],':k')
                  end
                  for KK=X2LO:X2HI
                      plot([X2FACE(KK) X2FACE(KK)],[X3FACE(X3LO) X3FACE(X3HI)],':k')
                  end
                  for IEDGE=CUT_EDGE_START+1:CUT_EDGE_END
                     NEDGE=MESHES(NM).CUT_EDGE(IEDGE).NEDGE;
                     for IEDGE2=1:NEDGE
                     n1  =MESHES(NM).CUT_EDGE(IEDGE).CEELEM(NOD1,IEDGE2);
                     n2  =MESHES(NM).CUT_EDGE(IEDGE).CEELEM(NOD2,IEDGE2);
                     XYZ1=MESHES(NM).CUT_EDGE(IEDGE).XYZVERT([X2AXIS X3AXIS],n1);
                     XYZ2=MESHES(NM).CUT_EDGE(IEDGE).XYZVERT([X2AXIS X3AXIS],n2);
                     plot([XYZ1(IAXIS) XYZ2(IAXIS)],[XYZ1(JAXIS) XYZ2(JAXIS)],'b','LineWidth',2);
                     end
                  end
                  pause
                  close(gcf)
               end

            end % I index
         end % J index
      end % K index

   end % X1AXIS_LOOP
   toc

   disp('Into GET_CARTCELL_CUTEDGES ...')
   tic
   % Now Define the INBOUNDARY cut-edge inside Cartesian cells:
   [ierr]=GET_CARTCELL_CUTEDGES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND);
   toc

   % 1. Cartesian GASPHASE cut-faces:
   % Loops for IAXIS, JAXIS, KAXIS faces: For FCVAR i,j,k, axis
   % - Define Cartesian Boundary Edges indexes.
   % - From ECVAR(i,j,k,IDCE,axis) figure out Entries in CUT_EDGE (GASPHASE segs).
   % - From FCVAR(i,j,k,IDCE,axis) figure out entries in CUT_EDGE (INBOUNDCF segs).
   % - Reorder Edges, figure out if there are disjoint areas present.
   % - Load into CUT_FACE <=> FCVAR(i,j,k,IDCF,axis).
   disp('Into GET_CARTFACE_CUTFACES ...')
   tic
   [ierr]=GET_CARTFACE_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,true);
   toc

   % 2. INBOUNDARY cut-faces:
   disp('Into GET_CARTCELL_CUTFACES ...')
   tic
   [ierr]=GET_CARTCELL_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,true);
   toc

   % Guard-cell Cartesian GASPHASE and INBOUNDARY cut-faces:
   disp('Into Guard cell GET_CARTFACE_CUTFACES, GET_CARTCELL_CUTFACES ...')
   tic
   [ierr]=GET_CARTFACE_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,false);
   [ierr]=GET_CARTCELL_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,false);
   toc

   % Finally: Definition of cut-cells:
   % Reset:
   CELLRT=zeros(IEND,JEND,KEND); % false
   MESHES(NM).N_SPCELL_CF = MESHES(NM).N_SPCELL;
   disp('Into GET_CARTCELL_CUTCELLS ...')
   tic
   [ierr]=GET_CARTCELL_CUTCELLS(NM);
   toc
end


% Loop over geometry:
SLEN_GEOM = 0.; AREA_GEOM = 0.; VOLUME_GEOM = 0.; XYZCEN_GEOM(IAXIS:KAXIS) = 0.;
for IG=1:N_GEOMETRY

   % Add length of wet surface edges:
   for IEDGE=1:GEOM(IG).N_EDGES
      SEG(NOD1:NOD2)  = GEOM(IG).EDGES(NOD1:NOD2,IEDGE);
      DV(IAXIS:KAXIS) = GEOM(IG).VERTS(MAX_DIM*(SEG(NOD2)-1)+1:MAX_DIM*SEG(NOD2)) - ...
                        GEOM(IG).VERTS(MAX_DIM*(SEG(NOD1)-1)+1:MAX_DIM*SEG(NOD1));
      SLEN = sqrt( DV(IAXIS)^2. + DV(JAXIS)^2. + DV(KAXIS)^2. );
      SLEN_GEOM = SLEN_GEOM + SLEN;
   end

   % Add to wet surface Areas:
   for IFACE=1:GEOM(IG).N_FACES
      if ( GEOM(IG).FACES_AREA(IFACE) < GEOMEPS )
         disp(['GEOM FACE=' num2str(IFACE) ', AREA=' ,num2str(GEOM(IG).FACES_AREA(IFACE))])
      end
   end
   AREA_GEOM = AREA_GEOM + GEOM(IG).GEOM_AREA;

   % Add to GEOMETRY volume:
   VOLUME_GEOM = VOLUME_GEOM + GEOM(IG).GEOM_VOLUME;

   % Add to XYZCEN for geometries:
   XYZCEN_GEOM(IAXIS:KAXIS)= XYZCEN_GEOM(IAXIS:KAXIS) + GEOM(IG).GEOM_VOLUME * GEOM(IG).GEOM_XYZCEN(IAXIS:KAXIS);
end
if(N_GEOMETRY > 0); XYZCEN_GEOM(IAXIS:KAXIS)=XYZCEN_GEOM(IAXIS:KAXIS)/VOLUME_GEOM; end

% Loop over meshes:
NCUTFACE_INB = 0;
CF_AREA_INB=0.;
CC_VOLUME_INB=0.;
GP_VOLUME=0.;
DM_VOLUME = 0.;
DM_XYZCEN(IAXIS:KAXIS) = 0.;
CCGP_XYZCEN(IAXIS:KAXIS) = 0.;
for NM=1:NMESHES

   ILO_CELL = MESHES(NM).ILO_CELL;  % First internal cell index. See notes.
   IHI_CELL = MESHES(NM).IHI_CELL;  % Last internal cell index.
   JLO_CELL = MESHES(NM).JLO_CELL;  % First internal cell index. See notes.
   JHI_CELL = MESHES(NM).JHI_CELL;  % Last internal cell index.
   KLO_CELL = MESHES(NM).KLO_CELL;  % First internal cell index. See notes.
   KHI_CELL = MESHES(NM).KHI_CELL;  % Last internal cell index.

   for ICF1 = 1:MESHES(NM).N_CUTFACE_MESH
     if (MESHES(NM).CUT_FACE(ICF1).STATUS == IBM_INBOUNDARY)
        NFACE = MESHES(NM).CUT_FACE(ICF1).NFACE;
        CF_AREA_INB = CF_AREA_INB + sum(MESHES(NM).CUT_FACE(ICF1).AREA(1:NFACE));
     end
   end

   for ICC1 = 1:MESHES(NM).N_CUTCELL_MESH
      NCELL = MESHES(NM).CUT_CELL(ICC1).NCELL;
      for ICC2=1:NCELL
         CCGP_XYZCEN(IAXIS:KAXIS) = CCGP_XYZCEN(IAXIS:KAXIS) + MESHES(NM).CUT_CELL(ICC1).VOLUME(ICC2) * ...
                                                               MESHES(NM).CUT_CELL(ICC1).XYZCEN(IAXIS:KAXIS,ICC2)';
         if ( MESHES(NM).CUT_CELL(ICC1).VOLUME(ICC2) < -GEOMEPS)
             disp(['Cut-cell=' num2str([ICC1,ICC2])  ', VOL=' num2str(MESHES(NM).CUT_CELL(ICC1).VOLUME(ICC2))])
         end
      end
      CC_VOLUME_INB = CC_VOLUME_INB + sum(MESHES(NM).CUT_CELL(ICC1).VOLUME(1:NCELL));
   end

   % Regular gasphase cells:
   for K=KLO_CELL:KHI_CELL
      for J=JLO_CELL:JHI_CELL
         for I=ILO_CELL:IHI_CELL

            if ( MESHES(NM).CCVAR(I,J,K,IBM_CGSC) == IBM_GASPHASE)
              SLEN = MESHES(NM).DX(I)*MESHES(NM).DY(J)*MESHES(NM).DZ(K);
              GP_VOLUME = GP_VOLUME + SLEN;
              CCGP_XYZCEN(IAXIS) = CCGP_XYZCEN(IAXIS) + SLEN * MESHES(NM).XC(I);
              CCGP_XYZCEN(JAXIS) = CCGP_XYZCEN(JAXIS) + SLEN * MESHES(NM).YC(J);
              CCGP_XYZCEN(KAXIS) = CCGP_XYZCEN(KAXIS) + SLEN * MESHES(NM).ZC(K);
            end
         end
      end
   end

   % Domain volume:
   SLEN = (MESHES(NM).XF-MESHES(NM).XS)*(MESHES(NM).YF-MESHES(NM).YS)*(MESHES(NM).ZF-MESHES(NM).ZS);
   DM_VOLUME = DM_VOLUME + SLEN;
   % Domain Centroid:
   DM_XYZCEN(IAXIS) = DM_XYZCEN(IAXIS) + SLEN * 0.5*(MESHES(NM).XF+MESHES(NM).XS);
   DM_XYZCEN(JAXIS) = DM_XYZCEN(JAXIS) + SLEN * 0.5*(MESHES(NM).YF+MESHES(NM).YS);
   DM_XYZCEN(KAXIS) = DM_XYZCEN(KAXIS) + SLEN * 0.5*(MESHES(NM).ZF+MESHES(NM).ZS);
end

disp(' Computational Geometry: Sanity tests for cut-cell region')
disp([' GEOM Surf  Area=' num2str(AREA_GEOM,'%11.4e') ', InBoundary Cut-faces Area=' num2str(CF_AREA_INB,'%11.4e') ...
', Relative Difference=' num2str((AREA_GEOM-CF_AREA_INB)/(AREA_GEOM+GEOMEPS),'%11.4e')])

CCGP_XYZCEN(IAXIS:KAXIS) = CCGP_XYZCEN(IAXIS:KAXIS) / (CC_VOLUME_INB+GP_VOLUME);
DM_XYZCEN(IAXIS:KAXIS)   = DM_XYZCEN(IAXIS:KAXIS)   / DM_VOLUME;

disp([' GEOM Gas Volume=' num2str(DM_VOLUME-VOLUME_GEOM,'%11.4e') ', Cut/Regl Gas cells Volume=' ...
      num2str(GP_VOLUME+CC_VOLUME_INB,'%11.4e') ', Relative Difference=' ...
      num2str(((DM_VOLUME-VOLUME_GEOM)-(GP_VOLUME+CC_VOLUME_INB))/(DM_VOLUME-VOLUME_GEOM),'%11.4e')])

disp([' GEOM Centroid               =' num2str(XYZCEN_GEOM(IAXIS:KAXIS),'%12.4e')])
disp([' DOMAIN-GEOM Centroid        =' ...
num2str((DM_XYZCEN(IAXIS:KAXIS)*DM_VOLUME - XYZCEN_GEOM(IAXIS:KAXIS)*VOLUME_GEOM)/(DM_VOLUME-VOLUME_GEOM),'%12.4e')])
disp([' Cut/Regl Gas cells Centroid =' num2str(CCGP_XYZCEN(IAXIS:KAXIS),'%12.4e')])
disp([' Centroid Relative Difference=' num2str(CCGP_XYZCEN(IAXIS:KAXIS)-...
(DM_XYZCEN(IAXIS:KAXIS)*DM_VOLUME - XYZCEN_GEOM(IAXIS:KAXIS)*VOLUME_GEOM)/(DM_VOLUME-VOLUME_GEOM),'%12.4e')])

ierr=0;

return
