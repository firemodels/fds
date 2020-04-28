function [ierr]=GET_CARTCELL_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,BNDINT_FLAG)


global XFACE DXFACE XCELL DXCELL YFACE DYFACE YCELL DYCELL ZFACE DZFACE ZCELL DZCELL
global ILO_CELL IHI_CELL JLO_CELL JHI_CELL KLO_CELL KHI_CELL
global ILO_FACE IHI_FACE JLO_FACE JHI_FACE KLO_FACE KHI_FACE
global IJK_COUNT IJK_COUNTF CCGUARD GEOMEPS FCELL
global NOD1 NOD2 NOD3 NOD4 IBM_MAX_WSTRIANG_SEG IAXIS JAXIS KAXIS
global IBM_FGSC IBM_IDCE IBM_CGSC IBM_CUTCFE IBM_SOLID IBM_GASPHASE
global IBM_UNDEFINED IBM_INBOUNDARY IBM_MAXVERTS_FACE IBM_MAXCEELEM_FACE
global GEOM MESHES IBM_IDCF IBM_GS LOW_IND HIGH_IND MAX_DIM
global BODINT_PLANE IBM_INBOUNDCF
global CELLRT

global X1FACE X2FACE X3FACE DX1FACE DX2FACE DX3FACE DX1CELL DX2CELL DX3CELL
global X2LO X2HI X3LO X3HI
global X1CELL X2CELL X3CELL X2LO_CELL X2HI_CELL X3LO_CELL X3HI_CELL


global basedir

global ISPCELL SPCELL_LIST

global N_GEOMETRY

plotflg = false;

ierr=1;

INDVERTBOD(1:3)  = [ 1, 2, 6 ];
INDVERTBOD2(1:3) = [ 2, 1, 6 ];

SEG_CELL=zeros(NOD2+IBM_MAX_WSTRIANG_SEG+2,20);
SEG_FLAG=zeros(1,IBM_MAXCEELEM_FACE);

% Define which cells are cut-cell, and which are solid:
if (BNDINT_FLAG)
   MESHES(NM).N_SPCELL=0;
   MESHES(NM).SPCELL_LIST=[0 0 0]';
   IJK_COUNT = zeros(IEND,JEND,KEND);
   IJK_COUNTF= zeros(IEND,JEND,KEND,MAX_DIM);
   ILO = ILO_CELL; IHI = IHI_CELL;
   JLO = JLO_CELL; JHI = JHI_CELL;
   KLO = KLO_CELL; KHI = KHI_CELL;
else
   ILO = ILO_CELL-CCGUARD; IHI = IHI_CELL+CCGUARD;
   JLO = JLO_CELL-CCGUARD; JHI = JHI_CELL+CCGUARD;
   KLO = KLO_CELL-CCGUARD; KHI = KHI_CELL+CCGUARD;
end

% Loop on Cartesian cells, define cut cells and solid cells ISSO:
for K=KLO:KHI
   for J=JLO:JHI
      for I=ILO:IHI

         if (IJK_COUNT(I,J,K)); continue; end

         % Face type of bounding Cartesian faces:
         FSID_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_FGSC,JAXIS);
         FSID_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_FGSC,JAXIS);
         FSID_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_FGSC,KAXIS);
         FSID_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_FGSC,KAXIS);

         % For this cell check if no Cartesian boundary faces are IBM_CUTCFE:
         % If outcell1 is true -> All regular faces for this face:
         OUTCELL1 = (FSID_XYZ(LOW_IND ,IAXIS) ~= IBM_CUTCFE) && ...
                    (FSID_XYZ(HIGH_IND,IAXIS) ~= IBM_CUTCFE) && ...
                    (FSID_XYZ(LOW_IND ,JAXIS) ~= IBM_CUTCFE) && ...
                    (FSID_XYZ(HIGH_IND,JAXIS) ~= IBM_CUTCFE) && ...
                    (FSID_XYZ(LOW_IND ,KAXIS) ~= IBM_CUTCFE) && ...
                    (FSID_XYZ(HIGH_IND,KAXIS) ~= IBM_CUTCFE);

         % Test for cell with INB edges:
         % If outcell2 is true -> no INB Edges associated with this cell:
         OUTCELL2 = 1; %(MESHES(NM).CCVAR(I,J,K,IBM_IDCE) <= 0);

         % Drop if outcell1 & outcell2
         if (OUTCELL1 && OUTCELL2)
            if ( (FSID_XYZ(LOW_IND ,IAXIS) == IBM_SOLID) && ...
                 (FSID_XYZ(HIGH_IND,IAXIS) == IBM_SOLID) && ...
                 (FSID_XYZ(LOW_IND ,JAXIS) == IBM_SOLID) && ...
                 (FSID_XYZ(HIGH_IND,JAXIS) == IBM_SOLID) && ...
                 (FSID_XYZ(LOW_IND ,KAXIS) == IBM_SOLID) && ...
                 (FSID_XYZ(HIGH_IND,KAXIS) == IBM_SOLID) );
               MESHES(NM).CCVAR(I,J,K,IBM_CGSC) = IBM_SOLID;
            end
            continue
         end

         MESHES(NM).CCVAR(I,J,K,IBM_CGSC) = IBM_CUTCFE;

      end
   end
end


% First add edges stemming from triangles laying on gridline planes:
% Dump triangle aligned segments as cut-cell cut-edges, on face cases:
% BNDINT_COND : if (BNDINT_FLAG)
   % Do Loop for different x1 planes:
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
          X1FACE = XFACE; DX1FACE = DXFACE;
          X2FACE = YFACE; DX2FACE = DYFACE;
          X3FACE = ZFACE; DX3FACE = DZFACE;
          % x2 cell center parameters:
          X2LO_CELL = JLO_CELL-CCGUARD; X2HI_CELL = JHI_CELL+CCGUARD;
          X2CELL = YCELL; DX2CELL = DYCELL;
          % x3 cell center parameters:
          X3LO_CELL = KLO_CELL-CCGUARD; X3HI_CELL = KHI_CELL+CCGUARD;
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
          X1FACE = YFACE; DX1FACE = DYFACE;
          X2FACE = ZFACE; DX2FACE = DZFACE;
          X3FACE = XFACE; DX3FACE = DXFACE;
          % x2 cell center parameters:
          X2LO_CELL = KLO_CELL-CCGUARD; X2HI_CELL = KHI_CELL+CCGUARD;
          X2CELL = ZCELL; DX2CELL = DZCELL;
          % x3 cell center parameters:
          X3LO_CELL = ILO_CELL-CCGUARD; X3HI_CELL = IHI_CELL+CCGUARD;
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
          X1FACE = ZFACE; DX1FACE = DZFACE;
          X2FACE = XFACE; DX2FACE = DXFACE;
          X3FACE = YFACE; DX3FACE = DYFACE;
          % x2 cell center parameters:
          X2LO_CELL = ILO_CELL-CCGUARD; X2HI_CELL = IHI_CELL+CCGUARD;
          X2CELL = XCELL; DX2CELL = DXCELL;
          % x3 cell center parameters:
          X3LO_CELL = JLO_CELL-CCGUARD; X3HI_CELL = JHI_CELL+CCGUARD;
          X3CELL = YCELL; DX3CELL = DYCELL;
      end

      % Loop Slices:
      for K=KLO:KHI
         for J=JLO:JHI
            for I=ILO:IHI

               IJK(IAXIS:KAXIS) = [ I, J, K ];

               % Plane:
               X1PLN = X1FACE(IJK(X1AXIS));

               % Get intersection of body on plane defined by X1PLN, normal to X1AXIS:
               DX2_MIN = min(DX2CELL(X2LO_CELL:X2HI_CELL));
               DX3_MIN = min(DX3CELL(X3LO_CELL:X3HI_CELL));
               TRI_ONPLANE_ONLY = true;
               RAYTRACE_X2_ONLY =false;
               [ierr,BODINT_PLANE]=GET_BODINT_PLANE(X1AXIS,X1PLN,IJK(X1AXIS),PLNORMAL,X2AXIS,...
                                                    X3AXIS,DX2_MIN,DX3_MIN,X2LO,X2HI,X3LO,X3HI,X2FACE,X3FACE,TRI_ONPLANE_ONLY,RAYTRACE_X2_ONLY);
               % Test that there is an intersection:
               if ((BODINT_PLANE.NTRIS) == 0); continue; end

               % Drop if node locations outside block plane area:
               if ((X2FACE(X2LO)-max(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))) > GEOMEPS); continue; end
               if ((min(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))-X2FACE(X2HI)) > GEOMEPS); continue; end
               if ((X3FACE(X3LO)-max(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))) > GEOMEPS); continue; end
               if ((min(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))-X3FACE(X3HI)) > GEOMEPS); continue; end

               % Allocate triangles variables:
               BODINT_PLANE.X1NVEC =zeros(1,BODINT_PLANE.NTRIS);
               BODINT_PLANE.AINV   =zeros(2,2,BODINT_PLANE.NTRIS);

               % Triangles inverses:
               for ITRI=1:BODINT_PLANE.NTRIS

                  TRIS(NOD1:NOD3) = BODINT_PLANE.TRIS(NOD1:NOD3,ITRI)';

                  % This is local IAXIS:JAXIS
                  XYEL(IAXIS:JAXIS,NOD1) = [ BODINT_PLANE.XYZ(X2AXIS,TRIS(NOD1)), ...
                                             BODINT_PLANE.XYZ(X3AXIS,TRIS(NOD1))  ]';
                  XYEL(IAXIS:JAXIS,NOD2) = [ BODINT_PLANE.XYZ(X2AXIS,TRIS(NOD2)), ...
                                             BODINT_PLANE.XYZ(X3AXIS,TRIS(NOD2))  ]';
                  XYEL(IAXIS:JAXIS,NOD3) = [ BODINT_PLANE.XYZ(X2AXIS,TRIS(NOD3)), ...
                                             BODINT_PLANE.XYZ(X3AXIS,TRIS(NOD3))  ]';

                  % Test that x1-x2-x3 obeys right hand rule:
                  VAL = (XYEL(IAXIS,NOD2)-XYEL(IAXIS,NOD1)) * (XYEL(JAXIS,NOD3)-XYEL(JAXIS,NOD1))- ...
                        (XYEL(JAXIS,NOD2)-XYEL(JAXIS,NOD1)) * (XYEL(IAXIS,NOD3)-XYEL(IAXIS,NOD1));

                  if (VAL >= 0.)
                     BODINT_PLANE.X1NVEC(ITRI) = 1;
                  else
                     BODINT_PLANE.X1NVEC(ITRI) =-1;
                  end

                  % Transformation Matrix for this triangle in x2x3 plane:
                  if (BODINT_PLANE.X1NVEC(ITRI) < 0.) % Rotate node 2 and 3 locations
                     DUMMY(IAXIS:JAXIS)     = XYEL(IAXIS:JAXIS,NOD2)';
                     XYEL(IAXIS:JAXIS,NOD2) = XYEL(IAXIS:JAXIS,NOD3);
                     XYEL(IAXIS:JAXIS,NOD3) = DUMMY(IAXIS:JAXIS)';
                  end

                  % Inverse of Master to physical triangle transform matrix:
                  A_COEF = XYEL(IAXIS,NOD1) - XYEL(IAXIS,NOD3);
                  B_COEF = XYEL(IAXIS,NOD2) - XYEL(IAXIS,NOD3);
                  C_COEF = XYEL(JAXIS,NOD1) - XYEL(JAXIS,NOD3);
                  D_COEF = XYEL(JAXIS,NOD2) - XYEL(JAXIS,NOD3);
                  DENOM  = A_COEF * D_COEF - B_COEF * C_COEF;
                  BODINT_PLANE.AINV(1,1,ITRI) =  D_COEF / DENOM;
                  BODINT_PLANE.AINV(2,1,ITRI) = -C_COEF / DENOM;
                  BODINT_PLANE.AINV(1,2,ITRI) = -B_COEF / DENOM;
                  BODINT_PLANE.AINV(2,2,ITRI) =  A_COEF / DENOM;

               end


               % There are triangles aligned with this x1pln:
               % Run by Face:
               % First solid Faces: x1 Faces, Check where they lay:
               for KK=X3LO_CELL:X3HI_CELL
                  for JJ=X2LO_CELL:X2HI_CELL

                     % Face indexes:
                     INDXI(IAXIS:KAXIS) = [ IJK(X1AXIS), JJ, KK ]; % Local x1,x2,x3
                     INDIF = INDXI(XIAXIS);
                     INDJF = INDXI(XJAXIS);
                     INDKF = INDXI(XKAXIS);

                     if(IJK_COUNTF(INDIF,INDJF,INDKF,X1AXIS)); continue; end
                     IJK_COUNTF(INDIF,INDJF,INDKF,X1AXIS) = true;

                     if (MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_FGSC,X1AXIS) ~= IBM_GASPHASE )

                        FVERT(IAXIS:JAXIS,NOD1) = [ X2FACE(JJ-FCELL  ), X3FACE(KK-FCELL  ) ]';
                        FVERT(IAXIS:JAXIS,NOD2) = [ X2FACE(JJ-FCELL+1), X3FACE(KK-FCELL  ) ]';
                        FVERT(IAXIS:JAXIS,NOD3) = [ X2FACE(JJ-FCELL+1), X3FACE(KK-FCELL+1) ]';
                        FVERT(IAXIS:JAXIS,NOD4) = [ X2FACE(JJ-FCELL  ), X3FACE(KK-FCELL+1) ]';

                        % Get triangle face intersection:
                        CEI = MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS);

                        % Triangle - face intersection vertices and edges:
                        [INB_FLG,FNVERT,XYVERT,FNEDGE,CEELEM,INDSEG]=GET_TRIANG_FACE_INT(X2AXIS,X3AXIS,FVERT,CEI,NM);

                        % XYvert to XYZvert:
                        if ( INB_FLG )
                           XYZVERTF = 0.;
                           XYZVERTF(X1AXIS,1:FNVERT) = X1PLN;
                           XYZVERTF(X2AXIS,1:FNVERT) = XYVERT(IAXIS,1:FNVERT);
                           XYZVERTF(X3AXIS,1:FNVERT) = XYVERT(JAXIS,1:FNVERT);

                           % Test for edges inside SOLID Region:
                           SOLID_EDGE(1:FNEDGE) = false;
                           for IEDGE=1:FNEDGE

                               % No body associated with segment. Might not be needed.
                               IG = INDSEG(4,IEDGE);
                               if( IG < 1); continue; end
                               SEG(NOD1:NOD2)  = CEELEM(NOD1:NOD2,IEDGE);
                               XP(IAXIS:KAXIS) = 0.5*(XYZVERTF(IAXIS:KAXIS,SEG(NOD1)) + ...
                                                      XYZVERTF(IAXIS:KAXIS,SEG(NOD2)));

                               % Direction NP:
                               NP(IAXIS:KAXIS) = 0.;
                               for I_NP=1:INDSEG(1,IEDGE)
                                   ITRI = INDSEG(1+I_NP,IEDGE);
                                   NP(IAXIS:KAXIS) = NP(IAXIS:KAXIS) + ...
                                   GEOM(IG).FACES_NORMAL(IAXIS:KAXIS,ITRI)';
                               end
                               NP = NP/norm(NP);
                               [val,XAXIS] = max(abs(NP(IAXIS:KAXIS)));

                               % Perturb XP in the average normal NP
                               % direction:
                               if(INDSEG(1,IEDGE) > 1)
                                   XP = XP + 10.*GEOMEPS*NP;
                               end

                               [SOLID_EDGE(IEDGE)]=GET_IS_SOLID_3D(XAXIS,XP,INDIF,INDJF,INDKF);

                           end

                           % Now drop SEGS with SOLID_EDGE(IEDGE)=true:
                           COUNT = 0;
                           for IEDGE=1:FNEDGE
                               if(SOLID_EDGE(IEDGE)); continue; end
                               COUNT=COUNT+1;
                               CEELEM(NOD1:NOD2,COUNT) = ...
                                   CEELEM(NOD1:NOD2,IEDGE);
                               INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,COUNT) = ...
                                   INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,IEDGE);
                           end
                           CEELEM(NOD1:NOD2,COUNT+1:FNEDGE)=IBM_UNDEFINED;
                           INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,COUNT+1:FNEDGE)=IBM_UNDEFINED;
                           FNEDGE = COUNT;

                           % Here ADD nodes and vertices to what is already
                           % there:
                           if (CEI == 0) % We need a new entry in CUT_EDGE
                              CEI      = MESHES(NM).N_CUTEDGE_MESH + 1;
                              MESHES(NM).N_CUTEDGE_MESH = CEI;
                              MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS) = CEI;
                              MESHES(NM).CUT_EDGE(CEI).NVERT  = FNVERT;
                              MESHES(NM).CUT_EDGE(CEI).NEDGE  = FNEDGE;
                              MESHES(NM).CUT_EDGE(CEI).NEDGE1 = 0;
                              MESHES(NM).CUT_EDGE(CEI).IJK(1:MAX_DIM+2) = ...
                                                   [ INDIF, INDJF, INDKF, X1AXIS, IBM_GS ];
                              MESHES(NM).CUT_EDGE(CEI).STATUS = IBM_INBOUNDCF;
                              MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,1:FNVERT) = ...
                                                      XYZVERTF(IAXIS:KAXIS,1:FNVERT);
                              MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,1:FNEDGE)    = ...
                                                       CEELEM(NOD1:NOD2,1:FNEDGE);
                              MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,1:FNEDGE) = ...
                                                       INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,1:FNEDGE);
                           else

                              NVERT_AUX=MESHES(NM).CUT_EDGE(CEI).NVERT;
                              NEDGE_OLD=MESHES(NM).CUT_EDGE(CEI).NEDGE;
                              for IVERT=1:FNVERT
                                 XYZV(IAXIS:KAXIS) = XYZVERTF(IAXIS:KAXIS,IVERT)';
                                 [NVERT_AUX,INOD]=INSERT_FACE_VERT(XYZV,NM,CEI,NVERT_AUX);
                                 for IEDGE=1:FNEDGE
                                    if (CEELEM(NOD1,IEDGE)==IVERT); CEELEM(NOD1,IEDGE)=INOD; end
                                    if (CEELEM(NOD2,IEDGE)==IVERT); CEELEM(NOD2,IEDGE)=INOD; end
                                 end
                              end


                              COUNT = NEDGE_OLD;
                              for IEDGE=1:FNEDGE

                                 % Here find if edge is already counted.
                                 FOUND=false;
                                 for IEOLD=1:NEDGE_OLD
                                 if (MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IEOLD) ~= CEELEM(NOD1,IEDGE)); continue; end
                                 if (MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IEOLD) ~= CEELEM(NOD2,IEDGE)); continue; end
                                 if (MESHES(NM).CUT_EDGE(CEI).INDSEG(4,IEOLD)    ~= INDSEG(4,IEDGE)); continue; end
                                 FOUND=true;
                                 end
                                 for IEOLD=1:NEDGE_OLD
                                 if (MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IEOLD) ~= CEELEM(NOD1,IEDGE)); continue; end
                                 if (MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IEOLD) ~= CEELEM(NOD2,IEDGE)); continue; end
                                 if (MESHES(NM).CUT_EDGE(CEI).INDSEG(4,IEOLD)    ~= INDSEG(4,IEDGE)); continue; end
                                 FOUND=true;
                                 end
                                 if (FOUND); continue; end

                                 COUNT=COUNT+1;
                                 MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,COUNT) = CEELEM(NOD1:NOD2,IEDGE);
                                 MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,COUNT)= ...
                                                          INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,IEDGE);

                              end



                              MESHES(NM).CUT_EDGE(CEI).NVERT = NVERT_AUX;
                              MESHES(NM).CUT_EDGE(CEI).NEDGE = COUNT;
                           end
                        end

                     end
                  end
               end

            end % I
         end % J
      end % K

   end
% end BNDINT_COND

% Second: Loop over cut-cells: For cut-cell i,j,k,lb
% - From cut-cell Cartesian faces, figure out INBOUNDCF segments (CUT_EDGE)
% and the wet surface triangles related to them.
% - From CCVAR(I,J,K,IBM_IDCE), firgure out INBOUNDCC segments in CUT_EDGE
% and triangles they belong to.
% - Working by triangle -> reorient segments using triangle normal outside
% of body (no disjoint areas are expected)
% - Load into CUT_FACE <=> CCVAR(I,J,K,IBM_IDCF).
if (BNDINT_FLAG)
   ILO = ILO_CELL; IHI = IHI_CELL;
   JLO = JLO_CELL; JHI = JHI_CELL;
   KLO = KLO_CELL; KHI = KHI_CELL;
else
   ILO = ILO_CELL-CCGUARD; IHI = IHI_CELL+CCGUARD;
   JLO = JLO_CELL-CCGUARD; JHI = JHI_CELL+CCGUARD;
   KLO = KLO_CELL-CCGUARD; KHI = KHI_CELL+CCGUARD;
end
% Loop on Cartesian cells, define cut cells and solid cells IBM_CGSC:
for K=KLO:KHI
   for J=JLO:JHI
      for I=ILO:IHI

         if ( MESHES(NM).CCVAR(I,J,K,IBM_CGSC) ~= IBM_CUTCFE ); continue; end

         if (CELLRT(I,J,K)); continue; end; % Special cell with bod-bod or self intersection.

         if (IJK_COUNT(I,J,K)); continue; end; IJK_COUNT(I,J,K)=true;

         % Face type of bounding Cartesian faces:
         FSID_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_FGSC,JAXIS);
         FSID_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_FGSC,JAXIS);
         FSID_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_FGSC,KAXIS);
         FSID_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_FGSC,KAXIS);

         % Start cut-cell INB cut-faces computation:
         % Loop local arrays to cell:
         NSEG      = 0;
         SEG_CELL  = IBM_UNDEFINED;

         NVERT     = 0;
         NFACE     = 0;
         XYZVERT   = 0.;

         % CUT_EDGE index of bounding Cartesian faces:
         CEIB_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_IDCE,IAXIS);
         CEIB_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_IDCE,IAXIS);
         CEIB_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_IDCE,JAXIS);
         CEIB_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_IDCE,JAXIS);
         CEIB_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_IDCE,KAXIS);
         CEIB_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_IDCE,KAXIS);

         % Cartesian Faces INBOUNDARY segments:
         for FAXIS=IAXIS:KAXIS
            for ILH=LOW_IND:HIGH_IND
               % By segment: Add Vertices/Segments to local arrays:
               CEI = CEIB_XYZ(ILH,FAXIS);
               if ( CEI > 0 ) % There are inboundary cut-edges
                  NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                  for IEDGE=1:NEDGE

                     SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE)';
                     STRI(1:IBM_MAX_WSTRIANG_SEG+2) = ...
                     MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,IEDGE)';

                     % x,y,z of node 1:
                     XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1))';
                     [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);
                     % x,y,z of node 2:
                     XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2))';
                     [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);

                     VEC(NOD1:NOD2) = [ INOD1, INOD2 ];
                     VEC(NOD2+1:NOD2+IBM_MAX_WSTRIANG_SEG+2) = STRI(1:IBM_MAX_WSTRIANG_SEG+2);
                     % Insertion ADD segment:
                     INLIST = false;
                     for IDUM = 1:NSEG
                        for IEQ1=1:3
                           EQUAL1 = SEG_CELL(INDVERTBOD(IEQ1),IDUM) == VEC(INDVERTBOD(IEQ1));
                           if (~EQUAL1); break; end
                        end
                        for IEQ2=1:3
                           EQUAL2 = SEG_CELL(INDVERTBOD(IEQ2),IDUM) == VEC(INDVERTBOD2(IEQ2));
                           if (~EQUAL2); break; end
                        end
                        if ( EQUAL1 || EQUAL2 )
                           if ( SEG_CELL(3,IDUM) > VEC(3) )
                              % for NOTHING:
                           elseif (SEG_CELL(3,IDUM) < VEC(3))
                              SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,IDUM) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
                           elseif (SEG_CELL(4,IDUM) ~= VEC(4))
                              SEG_CELL(3,IDUM) = SEG_CELL(3,IDUM) + 1;
                              SEG_CELL(5,IDUM) = VEC(4);
                           end
                           INLIST = true;
                           break
                        end
                     end
                     if (~INLIST)
                         NSEG = NSEG + 1;
                         SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,NSEG) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
                     end
                  end
               end
            end
         end

         % Cells INBOUNDARY segments:
         CEI = MESHES(NM).CCVAR(I,J,K,IBM_IDCE);
         if ( CEI > 0 ) % There are inboundary cut-edges
            NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
            for IEDGE=1:NEDGE

               SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE)';
               STRI(1:IBM_MAX_WSTRIANG_SEG+2) = ...
               MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,IEDGE)';

               % x,y,z of node 1:
               XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1))';
               [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);
               % x,y,z of node 2:
               XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2))';
               [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);

               if (INOD1 == INOD2); continue; end

               VEC(NOD1:NOD2) = [ INOD1, INOD2 ];
               VEC(NOD2+1:NOD2+IBM_MAX_WSTRIANG_SEG+2) = STRI(1:IBM_MAX_WSTRIANG_SEG+2);
               % Insertion ADD segment:
               INLIST = false;
               for IDUM = 1:NSEG
                  for IEQ1=1:3
                     EQUAL1 = SEG_CELL(INDVERTBOD(IEQ1),IDUM) == VEC(INDVERTBOD(IEQ1));
                     if (~EQUAL1); break; end
                  end
                  for IEQ2=1:3
                     EQUAL2 = SEG_CELL(INDVERTBOD(IEQ2),IDUM) == VEC(INDVERTBOD2(IEQ2));
                     if (~EQUAL2); break; end
                  end
                  if ( EQUAL1 || EQUAL2 )
                     if ( SEG_CELL(3,IDUM) > VEC(3) )
                        % for NOTHING:
                     elseif (SEG_CELL(3,IDUM) < VEC(3))
                        SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,IDUM) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
                     elseif (SEG_CELL(4,IDUM) ~= VEC(4))
                        SEG_CELL(3,IDUM) = SEG_CELL(3,IDUM) + 1;
                        SEG_CELL(5,IDUM) = VEC(4);
                     end
                     INLIST = true;
                     break
                  end
               end
               if (~INLIST)
                   NSEG = NSEG + 1;
                   SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,NSEG) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
               end
            end
         end



         % Drop segments that are unconnected:
         VERT_SEGS = zeros(1,NVERT);
         for IDUM = 1:NSEG
         VERT_SEGS(SEG_CELL(NOD1,IDUM)) = VERT_SEGS(SEG_CELL(NOD1,IDUM))+1;
         VERT_SEGS(SEG_CELL(NOD2,IDUM)) = VERT_SEGS(SEG_CELL(NOD2,IDUM))+1;
         end
         SEG_CELL2 = SEG_CELL;
         COUNT = 0;
         for IDUM = 1:NSEG
         if(VERT_SEGS(SEG_CELL2(NOD1,IDUM)) > 1 && VERT_SEGS(SEG_CELL2(NOD2,IDUM)) > 1)
             COUNT = COUNT + 1;
             SEG_CELL(:,COUNT) = SEG_CELL2(:,IDUM);
             continue
         end
         end
         NSEG = COUNT;


%          if(I==25 && J==22 && K==8)
%              disp('1 Found cell 25,22,8')
%              SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,1:NSEG)
%              XYZVERT(:,1:NVERT)'
%              figure
%              hold on
%              a=0.001;
%              colr = ['-or'; '-ob'; '-om'];
%              for ISEG=1:NSEG
%                  XYZSEG(:,NOD1)=XYZVERT(:,SEG_CELL(NOD1,ISEG));
%                  XYZSEG(:,NOD2)=XYZVERT(:,SEG_CELL(NOD2,ISEG));
%                  IG         =max(1,SEG_CELL(6,ISEG));
%                  plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),colr(IG,1:3),'LineWidth',2)
%                  text(0.5*(sum(XYZSEG(IAXIS,:)))+a,0.5*(sum(XYZSEG(JAXIS,:)))+a,...
%                       0.5*(sum(XYZSEG(KAXIS,:)))+a,num2str(SEG_CELL(4,ISEG)),'FontSize',10)
%              end
%              for IVERT=1:NVERT
%                  text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                       XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
%
%              end
%              xlabel('X')
%              ylabel('Y')
%              zlabel('Z')
%
%              P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
%              P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
%              P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
%              P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
%              P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
%              P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
%              P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
%              P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];
%
%              FC(1,:) = [ 1 4 8 5 ];
%              FC(2,:) = [ 2 3 7 6 ];
%              FC(3,:) = [ 1 5 6 2 ];
%              FC(4,:) = [ 4 8 7 3 ];
%              FC(5,:) = [ 1 2 3 4 ];
%              FC(6,:) = [ 5 6 7 8 ];
%
%              for IFC=1:6
%                  [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
%                  set(hp,'FaceAlpha',0.3)
%              end
%
%              %axis equal
%              box on
%              view([45 45])
%              pause
%          end


         % Now obtain body-triangle combinations present:
         BOD_TRI = IBM_UNDEFINED;
         NBODTRI = 0;
         for ISEG=1:NSEG

            % First triangle location (Assume one body and at
            % most two triangs per seg).
            INLIST = false;
            for IBODTRI=1:NBODTRI
               if ( (BOD_TRI(1,IBODTRI) == SEG_CELL(6,ISEG)) && ...
                    (BOD_TRI(2,IBODTRI) == SEG_CELL(4,ISEG)) )
                  % Body/triang already on list.
                  INLIST = true;
                  continue
               end
            end
            if (~INLIST)
               % Add first triang to list:
               NBODTRI = NBODTRI + 1;
               BOD_TRI(1:2,NBODTRI) = SEG_CELL( [ 6, 4 ] , ISEG);
            end

            % No second triangle associated:
            if ( SEG_CELL(3,ISEG) < 2 ); continue; end

            % Second triangle location
            INLIST = false;
            for IBODTRI=1:NBODTRI
               if ( (BOD_TRI(1,IBODTRI) == SEG_CELL(6,ISEG)) && ...
                    (BOD_TRI(2,IBODTRI) == SEG_CELL(5,ISEG)) )
                  % Body/triang already on list.
                  INLIST = true;
                  continue
               end
            end
            if (~INLIST)
               % Add first triang to list:
               NBODTRI = NBODTRI + 1;
               BOD_TRI(1:2,NBODTRI) = SEG_CELL( [ 6, 5 ] , ISEG);
            end
         end % ISEG.

         % Do Test for cycling when all body-triangle combinations produce two or less segments:
         SEG_FLAG(1)=true;
         for ICF=1:NBODTRI
            IBOD = BOD_TRI(1,ICF);
            ITRI = BOD_TRI(2,ICF);
            NSEG_FACE = 0;
            for ISEG=1:NSEG
               if ((SEG_CELL(6,ISEG) == IBOD) && ...
                  ((SEG_CELL(4,ISEG) == ITRI) || (SEG_CELL(5,ISEG) == ITRI)) )
                  NSEG_FACE = NSEG_FACE + 1;
               end
            end
            % If only one or two seg => continue:
            if ( NSEG_FACE <= 2 ); continue; end
            SEG_FLAG(1)=false;
            break
         end
         if (SEG_FLAG(1)); continue; end % CYCLES I,J,K loop.

         % This is a cut-face, allocate space:
         NCUTFACE = MESHES(NM).N_CUTFACE_MESH + MESHES(NM).N_GCCUTFACE_MESH + 1;
         if (BNDINT_FLAG)
            MESHES(NM).N_CUTFACE_MESH = NCUTFACE;
         else
            MESHES(NM).N_GCCUTFACE_MESH = MESHES(NM).N_GCCUTFACE_MESH + 1;
         end
         MESHES(NM).CCVAR(I,J,K,IBM_IDCF) = NCUTFACE;
         MESHES(NM).CUT_FACE(NCUTFACE).NVERT  = NVERT;
         MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = 0;
         MESHES(NM).CUT_FACE(NCUTFACE).IJK(1:MAX_DIM+1) = [ I, J, K, 0 ]; % No axis = 0
         MESHES(NM).CUT_FACE(NCUTFACE).STATUS = IBM_INBOUNDARY;
         MESHES(NM).CUT_FACE(NCUTFACE).XYZVERT(IAXIS:KAXIS,1:NVERT) = XYZVERT(IAXIS:KAXIS,1:NVERT);

         % Running by body-triangle combination, define list of
         % segments that belong to each pair.
         for ICF=1:NBODTRI

            IBOD = BOD_TRI(1,ICF);
            ITRI = BOD_TRI(2,ICF);

            SEG_FACE  = IBM_UNDEFINED;
            NSEG_FACE = 0;
            for ISEG=1:NSEG
               if ((SEG_CELL(6,ISEG) == IBOD) && ...
                  ((SEG_CELL(4,ISEG) == ITRI) || (SEG_CELL(5,ISEG) == ITRI)) )
                  NSEG_FACE = NSEG_FACE + 1;
                  SEG_FACE(NOD1:NOD2,NSEG_FACE) = SEG_CELL(NOD1:NOD2,ISEG);
               end
            end

            % If only one or two seg => continue:
            if ( NSEG_FACE <= 2 ); continue; end

            % Now build sequential list of segments:
            SEG_FACE2 = IBM_UNDEFINED*ones(2,NSEG_FACE); %[nod1 nod2]
            SEG_FLAG  = ones(1,NSEG_FACE);
            ISEG_FACE = 1;
            COUNTR    = 1;
            CTSTART   = COUNTR;
            SEG_FACE2(NOD1:NOD2,COUNTR) = SEG_FACE(NOD1:NOD2,ISEG_FACE);
            SEG_FLAG(ISEG_FACE) = false;
            NSEG_LEFT = NSEG_FACE - 1;
            CTR = 0;
            % Infinite Loop:
            CYCLE_CELL=false;
            while(1)
               for ISEG_FACE=1:NSEG_FACE

                  if (SEG_FLAG(ISEG_FACE)) % This seg hasn't been added to seg_face2
                     % Test for common node:
                     if ( SEG_FACE2(NOD2,COUNTR) == SEG_FACE(NOD1,ISEG_FACE) )
                        COUNTR = COUNTR + 1;
                        SEG_FACE2(NOD1:NOD2,COUNTR) = SEG_FACE(NOD1:NOD2,ISEG_FACE);
                        SEG_FLAG(ISEG_FACE) = false;
                        NSEG_LEFT = NSEG_LEFT - 1;
                        break
                     elseif ( SEG_FACE2(NOD2,COUNTR) == SEG_FACE(NOD2,ISEG_FACE) )

                        if ( SEG_FACE2(NOD1,COUNTR) == SEG_FACE(NOD1,ISEG_FACE) )
                           disp('Building INBOUND faces, repeated index.')
                        end
                        COUNTR = COUNTR + 1;
                        SEG_FACE2(NOD1:NOD2,COUNTR) = SEG_FACE( [ NOD2, NOD1 ] ,ISEG_FACE);
                        SEG_FLAG(ISEG_FACE) = false;
                        NSEG_LEFT = NSEG_LEFT - 1;
                        break
                     end
                  end
               end
               % Break loop:
               if ( NSEG_LEFT == 0 ); break; end
               CTR = CTR + 1;

               % Plot cell and cut-faces if there is no convergence:
               if ( CTR > NSEG_FACE^3 )

%                       disp(['Error GET_CARTCELL_CUTFACES: ctr > nseg_face^3 ,' num2str(BNDINT_FLAG) ',' ...
%                            num2str(I) ',' num2str(J) ',' num2str(K) ',' num2str(NCUTFACE) ',' ...
%                            num2str(MESHES(NM).CUT_FACE(NCUTFACE).NFACE)])
                      CYCLE_CELL=true;
                      MESHES(NM).N_SPCELL = MESHES(NM).N_SPCELL + 1;
                      MESHES(NM).SPCELL_LIST(IAXIS:KAXIS,MESHES(NM).N_SPCELL) = [I J K]';



%                         file=[basedir '/Cartcell_cutfaces_' num2str(I) '_' num2str(J) '_' num2str(K) '.dat'];
%                         %disp(['Writing file : ' file])
%                         [fid]=fopen(file,'w');
%
%                         % Info pertaining to the Cartesian Cell:
%                         fprintf(fid,'I,J,K:\n');
%                         fprintf(fid,'%d  %d  %d  %0.6g\n',I,J,K,GEOMEPS);
%                         fprintf(fid,'XC(I),DX(I),YC(J),DY(J),ZC(K),DZ(K)\n');
%                         fprintf(fid,'%12.8f %12.8f\n',[XCELL(I) DXCELL(I)]); % MESHES(NM).XC(I),MESHES(NM).DX(I)
%                         fprintf(fid,'%12.8f %12.8f\n',[YCELL(J) DYCELL(J)]); % MESHES(NM).YC(J),MESHES(NM).DY(J)
%                         fprintf(fid,'%12.8f %12.8f\n',[ZCELL(K) DZCELL(K)]); % MESHES(NM).ZC(K),MESHES(NM).DZ(K)
%                         fprintf(fid,'NVERT,NSEG,NSEG_FACE,COUNTR,NSEG_LEFT:\n');
%                         fprintf(fid,'%d %d %d %d %d\n',[NVERT,NSEG,NSEG_FACE,COUNTR,NSEG_LEFT]);
%                         fprintf(fid,'XYZVERT(IAXIS:KAXIS,1:NVERT):\n');
%                         for IDUM=1:NVERT
%                            fprintf(fid,'%d %12.8f %12.8f %12.8f\n',[IDUM XYZVERT(IAXIS:KAXIS,IDUM)']);
%                         end
%                         fprintf(fid,'SEG_CELL(NOD1:NOD2,1:NSEG),SEG_CELL(3:6,1:NSEG):\n');
%                         for IDUM=1:NSEG
%                            fprintf(fid,'%d %d %d %d %d %d %d\n',[IDUM SEG_CELL(NOD1:NOD2,IDUM)' SEG_CELL(3:6,IDUM)']);
%                         end
%                         fprintf(fid,'SEG_FACE(NOD1:NOD2,1:NSEG_FACE):\n');
%                         for IDUM=1:NSEG_FACE
%                            fprintf(fid,'%d %d %d\n',[IDUM SEG_FACE(NOD1:NOD2,IDUM)']);
%                         end
%                         fprintf(fid,'SEG_FACE2(NOD1:NOD21:COUNTR):\n');
%                         for IDUM=1:COUNTR
%                            fprintf(fid,'%d %d %d\n',[IDUM SEG_FACE2(NOD1:NOD2,IDUM)']);
%                         end
%                         fprintf(fid,'ICF,BOD_TRI:\n');
%                         fprintf(fid,'%d %d\n',[ICF NBODTRI]);
%                         for IDUM=1:NBODTRI
%                            fprintf(fid,'%d %d\n',BOD_TRI(1:2,IDUM)');
%                         end
%                         fclose(fid);
%                         %disp(['File written.'])
                    break

               end

            end
            if(CYCLE_CELL); break; end

            if ( COUNTR ~= NSEG_FACE)
                disp('Building INBOUND faces: ~isequal(countr,nseg)')
            end

            % Using triangles normal, reorder nodes as in right hand rule.
            NORMTRI(IAXIS:KAXIS) = GEOM(IBOD).FACES_NORMAL(IAXIS:KAXIS,ITRI);

            % First test if INB face is on Cartesian face and pointing
            % outside of Cartesian cell. If so drop:
            XYZ(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,SEG_FACE2(1,1))';
            % IAXIS:
            if ( (abs(NORMTRI(IAXIS)+1.) < GEOMEPS) && ...
                 (abs(XFACE(I-FCELL  )-XYZ(IAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(IAXIS)-1.) < GEOMEPS) && ...
                 (abs(XFACE(I-FCELL+1)-XYZ(IAXIS)) < GEOMEPS) ); continue; end % High Face
            % JAXIS:
            if ( (abs(NORMTRI(JAXIS)+1.) < GEOMEPS) && ...
                 (abs(YFACE(J-FCELL  )-XYZ(JAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(JAXIS)-1.) < GEOMEPS) && ...
                 (abs(YFACE(J-FCELL+1)-XYZ(JAXIS)) < GEOMEPS) ); continue; end % High Face
            % KAXIS:
            if ( (abs(NORMTRI(KAXIS)+1.) < GEOMEPS) && ...
                 (abs(ZFACE(K-FCELL  )-XYZ(KAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(KAXIS)-1.) < GEOMEPS) && ...
                 (abs(ZFACE(K-FCELL+1)-XYZ(KAXIS)) < GEOMEPS) ); continue; end % High Face

            % Face Vertices average location:
            XCEN(IAXIS:KAXIS) = 0.;
            for ISEG_FACE=1:NSEG_FACE
                XCEN(IAXIS:KAXIS) = XCEN(IAXIS:KAXIS) + XYZVERT(IAXIS:KAXIS,SEG_FACE2(NOD1,ISEG_FACE))';
            end
            XCEN(IAXIS:KAXIS) = XCEN(IAXIS:KAXIS) / NSEG_FACE;

            ISEG_FACE = 1;
            VC1(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,SEG_FACE2(NOD1,ISEG_FACE  ))' - XCEN(IAXIS:KAXIS);
            V12(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,SEG_FACE2(NOD1,ISEG_FACE+1))' - ...
                               XYZVERT(IAXIS:KAXIS,SEG_FACE2(NOD1,ISEG_FACE  ))';

            CROSSV(IAXIS) = VC1(JAXIS)*V12(KAXIS) - VC1(KAXIS)*V12(JAXIS);
            CROSSV(JAXIS) = VC1(KAXIS)*V12(IAXIS) - VC1(IAXIS)*V12(KAXIS);
            CROSSV(KAXIS) = VC1(IAXIS)*V12(JAXIS) - VC1(JAXIS)*V12(IAXIS);

            RH_ORIENTED = ( NORMTRI(IAXIS)*CROSSV(IAXIS) + ...
                            NORMTRI(JAXIS)*CROSSV(JAXIS) + ...
                            NORMTRI(KAXIS)*CROSSV(KAXIS) ) > 0.;

            NP  = NSEG_FACE;
            NCF = MESHES(NM).CUT_FACE(NCUTFACE).NFACE + 1;
%             NVSIZE=SIZE(MESHES(NM).CUT_FACE(NCUTFACE).CFELEM,DIM=1)
%             if (NP+1 > NVSIZE)
%                % REALLOCATE CFELEM:
%                ALLOCATE(CFELEM(1:NVSIZE,1:NBODTRI))
%                CFELEM(1:NVSIZE,1:NBODTRI) = MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1:NVSIZE,1:NBODTRI)
%                DEALLOCATE(MESHES(NM).CUT_FACE(NCUTFACE).CFELEM)
%                ALLOCATE(MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1:NP+1+DELTA_VERT,1:NBODTRI))
%                MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(:,:) = IBM_UNDEFINED
%                MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1:NVSIZE,1:NBODTRI) = CFELEM(1:NVSIZE,1:NBODTRI)
%                DEALLOCATE(CFELEM)
%             end
            MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1,NCF) = NP;
            if (RH_ORIENTED)
                for IDUM=1:NP
                   MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(IDUM+1,NCF) = SEG_FACE2(NOD1,IDUM);
                end
            else
                for IDUM=1:NP
                   MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(IDUM+1,NCF) = SEG_FACE2(NOD1,NP+1-IDUM);
                end
            end
            MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = NCF;

            % Compute Sections area and centroid:
            AREA                = 0.;
            ACEN(IAXIS:KAXIS)   = 0.;
            INXAREA             = 0.;
            SQAREA(IAXIS:KAXIS) = 0.;
            for ISEG_FACE=1:NSEG_FACE-1

               IDUM = MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1+ISEG_FACE,NCF);
               X1(IAXIS:KAXIS)  = XYZVERT(IAXIS:KAXIS,IDUM);
               IDUM = MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(2+ISEG_FACE,NCF);
               X2(IAXIS:KAXIS)  = XYZVERT(IAXIS:KAXIS,IDUM);
               VC1(IAXIS:KAXIS) = X1(IAXIS:KAXIS) - XCEN(IAXIS:KAXIS);
               V12(IAXIS:KAXIS) = X2(IAXIS:KAXIS) - X1(IAXIS:KAXIS);
               XCENI(IAXIS:KAXIS) = (XCEN(IAXIS:KAXIS) + X1(IAXIS:KAXIS) + X2(IAXIS:KAXIS)) / 3.;

               CROSSV(IAXIS) = VC1(JAXIS)*V12(KAXIS) - VC1(KAXIS)*V12(JAXIS);
               CROSSV(JAXIS) = VC1(KAXIS)*V12(IAXIS) - VC1(IAXIS)*V12(KAXIS);
               CROSSV(KAXIS) = VC1(IAXIS)*V12(JAXIS) - VC1(JAXIS)*V12(IAXIS);

               AREAI = 0.5 * sqrt( CROSSV(IAXIS)^2. + CROSSV(JAXIS)^2. + CROSSV(KAXIS)^2. );
               AREA  = AREA + AREAI;
               ACEN(IAXIS:KAXIS)  = ACEN(IAXIS:KAXIS) + AREAI * XCENI(IAXIS:KAXIS);
               % volume computation variables:
               XC1(IAXIS:KAXIS) = 0.5*(XCEN(IAXIS:KAXIS) + X1(IAXIS:KAXIS));
               XC2(IAXIS:KAXIS) = 0.5*(XCEN(IAXIS:KAXIS) + X2(IAXIS:KAXIS));
               X12(IAXIS:KAXIS) = 0.5*(  X1(IAXIS:KAXIS) + X2(IAXIS:KAXIS));
               % dot(i,nc) int(x)dA
               INXAREA = INXAREA + NORMTRI(IAXIS)*XCENI(IAXIS)*AREAI; % Single Gauss pt integration.
               % dot(i,nc) int(x^2)dA, dot(j,nc) int(y^2)dA, dot(k,nc) int(z^2)dA
               for IX=IAXIS:KAXIS
                  INT2 = (XC1(IX)^2. + XC2(IX)^2. + X12(IX)^2.) / 3.;
                  SQAREA(IX) = SQAREA(IX) + NORMTRI(IX)*INT2*AREAI;  % Midpoint rule.
               end
            end
            % Final seg:
            IDUM = MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1+NSEG_FACE,NCF);
            X1(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,IDUM);
            IDUM = MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1+1        ,NCF);
            X2(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,IDUM);

            VC1(IAXIS:KAXIS) = X1(IAXIS:KAXIS) - XCEN(IAXIS:KAXIS);
            V12(IAXIS:KAXIS) = X2(IAXIS:KAXIS) - X1(IAXIS:KAXIS);
            XCENI(IAXIS:KAXIS) = (XCEN(IAXIS:KAXIS) + X1(IAXIS:KAXIS) + X2(IAXIS:KAXIS)) / 3.;

            CROSSV(IAXIS) = VC1(JAXIS)*V12(KAXIS) - VC1(KAXIS)*V12(JAXIS);
            CROSSV(JAXIS) = VC1(KAXIS)*V12(IAXIS) - VC1(IAXIS)*V12(KAXIS);
            CROSSV(KAXIS) = VC1(IAXIS)*V12(JAXIS) - VC1(JAXIS)*V12(IAXIS);

            AREAI = 0.5 * sqrt( CROSSV(IAXIS)^2. + CROSSV(JAXIS)^2. + CROSSV(KAXIS)^2. );
            AREA  = AREA + AREAI;
            ACEN(IAXIS:KAXIS)  = (ACEN(IAXIS:KAXIS) + AREAI * XCENI(IAXIS:KAXIS))/AREA;
            % volume computation variables:
            XC1(IAXIS:KAXIS) = 0.5*(XCEN(IAXIS:KAXIS) + X1(IAXIS:KAXIS));
            XC2(IAXIS:KAXIS) = 0.5*(XCEN(IAXIS:KAXIS) + X2(IAXIS:KAXIS));
            X12(IAXIS:KAXIS) = 0.5*(  X1(IAXIS:KAXIS) + X2(IAXIS:KAXIS));
            % dot(i,nc) int(x)dA
            INXAREA = INXAREA + NORMTRI(IAXIS)*XCENI(IAXIS)*AREAI; % Single Gauss pt integration.
            % dot(i,nc) int(x^2)dA, dot(j,nc) int(y^2)dA, dot(k,nc) int(z^2)dA
            for IX=IAXIS:KAXIS
               INT2 = (XC1(IX)^2. + XC2(IX)^2. + X12(IX)^2.) / 3.;
               SQAREA(IX) = SQAREA(IX) + NORMTRI(IX)*INT2*AREAI;  % Midpoint rule.
            end

            MESHES(NM).CUT_FACE(NCUTFACE).AREA(NCF) = AREA;
            MESHES(NM).CUT_FACE(NCUTFACE).XYZCEN(IAXIS:KAXIS,NCF) = ACEN(IAXIS:KAXIS);

            % Fields for cut-cell volume/centroid computation:
            % dot(i,nc)*int(x)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).INXAREA(NCF)   = INXAREA;
            % dot(i,nc)*int(x^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).INXSQAREA(NCF) = SQAREA(IAXIS);
            % dot(j,nc)*int(y^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).JNYSQAREA(NCF) = SQAREA(JAXIS);
            % dot(k,nc)*int(z^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).KNZSQAREA(NCF) = SQAREA(KAXIS);

            % Define Body-triangle reference:
            MESHES(NM).CUT_FACE(NCUTFACE).BODTRI(1:2,NCF)= [ IBOD, ITRI ];

            % Assign surf-index: Depending on GEOMETRY:
            MESHES(NM).CUT_FACE(NCUTFACE).SURF_INDEX(NCF) = GEOM(IBOD).SURFS(ITRI);

         end

         if(CYCLE_CELL)

             CELLRT(I,J,K) = 1;
             IJK_COUNT(I,J,K)=false;

             MESHES(NM).CCVAR(I,J,K,IBM_IDCF) = IBM_UNDEFINED;
             MESHES(NM).CUT_FACE(NCUTFACE).NVERT  = 0;
             MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = 0;
             MESHES(NM).CUT_FACE(NCUTFACE).IJK(1:MAX_DIM+1) = 0; % No axis = 0
             MESHES(NM).CUT_FACE(NCUTFACE).STATUS = IBM_UNDEFINED;
             MESHES(NM).CUT_FACE(NCUTFACE).XYZVERT=0;
             MESHES(NM).CUT_FACE(NCUTFACE).CEELEM=0;

             % This is a cut-face, allocate space:
             NCUTFACE = NCUTFACE-1; %MESHES(NM).N_CUTFACE_MESH + MESHES(NM).N_GCCUTFACE_MESH + 1;
             if (BNDINT_FLAG)
                MESHES(NM).N_CUTFACE_MESH = NCUTFACE;
             else
                MESHES(NM).N_GCCUTFACE_MESH = MESHES(NM).N_GCCUTFACE_MESH - 1;
             end

         end
      end % I
   end % J
end % K

% Now process special cells of type CELLRT=T:
% Loop on Cartesian cells, define cut cells and solid cells IBM_CGSC:
for K=KLO:KHI
   for J=JLO:JHI
      for I=ILO:IHI

         if ( MESHES(NM).CCVAR(I,J,K,IBM_CGSC) ~= IBM_CUTCFE ); continue; end
         if (~CELLRT(I,J,K)); continue; end; % Special cell with bod-bod or self intersection.
         if (IJK_COUNT(I,J,K)); continue; end; IJK_COUNT(I,J,K)=true;

         % Start cut-cell INB cut-faces computation:
         % Loop local arrays to cell:
         NSEG      = 0;
         SEG_CELL  = IBM_UNDEFINED;

         NVERT     = 0;
         NFACE     = 0;
         XYZVERT   = 0.;

         % CUT_EDGE index of bounding Cartesian faces:
         CEIB_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_IDCE,IAXIS);
         CEIB_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_IDCE,IAXIS);
         CEIB_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_IDCE,JAXIS);
         CEIB_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_IDCE,JAXIS);
         CEIB_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_IDCE,KAXIS);
         CEIB_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_IDCE,KAXIS);

         % Cartesian Faces INBOUNDARY segments:
         for FAXIS=IAXIS:KAXIS
            for ILH=LOW_IND:HIGH_IND
               % By segment: Add Vertices/Segments to local arrays:
               CEI = CEIB_XYZ(ILH,FAXIS);
               if ( CEI > 0 ) % There are inboundary cut-edges
                  NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE1;
                  for IEDGE=1:NEDGE

                     SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE)';
                     STRI(1:IBM_MAX_WSTRIANG_SEG+2) = ...
                     MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,IEDGE)';

                     % x,y,z of node 1:
                     XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1))';
                     [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);
                     % x,y,z of node 2:
                     XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2))';
                     [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZ,NVERT,XYZVERT);

                     VEC(NOD1:NOD2) = (HIGH_IND-ILH)*[ INOD1, INOD2 ] + (ILH-LOW_IND)*[ INOD2, INOD1 ];
                     VEC(NOD2+1:NOD2+IBM_MAX_WSTRIANG_SEG+2) = STRI(1:IBM_MAX_WSTRIANG_SEG+2);

                     % Insertion ADD segment:
                     INLIST = false;
                     for IDUM = 1:NSEG
                        for IEQ1=1:3
                           EQUAL1 = SEG_CELL(INDVERTBOD(IEQ1),IDUM) == VEC(INDVERTBOD(IEQ1));
                           if (~EQUAL1); break; end
                        end
                        for IEQ2=1:3
                           EQUAL2 = SEG_CELL(INDVERTBOD(IEQ2),IDUM) == VEC(INDVERTBOD2(IEQ2));
                           if (~EQUAL2); break; end
                        end
                        if ( EQUAL1 || EQUAL2 )
                           if ( SEG_CELL(3,IDUM) > VEC(3) )
                              % for NOTHING:
                           elseif (SEG_CELL(3,IDUM) < VEC(3))
                              SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,IDUM) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
                           elseif (SEG_CELL(4,IDUM) ~= VEC(4))
                              SEG_CELL(3,IDUM) = SEG_CELL(3,IDUM) + 1;
                              SEG_CELL(5,IDUM) = VEC(4);
                           end
                           INLIST = true;
                           break
                        end
                     end
                     if (~INLIST)
                         NSEG = NSEG + 1;
                         SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,NSEG) = VEC(1:NOD2+IBM_MAX_WSTRIANG_SEG+2);
                         SEG_POS(NSEG) = (2*ILH-3)*FAXIS;
                     end

                  end
               end
            end
         end

         % Drop segments that are unconnected:
         VERT_SEGS = zeros(1,NVERT);
         for IDUM = 1:NSEG
         VERT_SEGS(SEG_CELL(NOD1,IDUM)) = VERT_SEGS(SEG_CELL(NOD1,IDUM))+1;
         VERT_SEGS(SEG_CELL(NOD2,IDUM)) = VERT_SEGS(SEG_CELL(NOD2,IDUM))+1;
         end
         SEG_CELL2 = SEG_CELL;
         COUNT = 0;
         for IDUM = 1:NSEG
         if(VERT_SEGS(SEG_CELL2(NOD1,IDUM)) > 1 && VERT_SEGS(SEG_CELL2(NOD2,IDUM)) > 1)
             COUNT = COUNT + 1;
             SEG_CELL(:,COUNT) = SEG_CELL2(:,IDUM);
             continue
         end
         end
         NSEG = COUNT;

%          if(I==36 && J==17 && K==24)
%              disp('2 Found cell 36,17,24')
%              SEG_CELL(1:NOD2+IBM_MAX_WSTRIANG_SEG+2,1:NSEG)
%              XYZVERT(:,1:NVERT)'
%              figure
%              hold on
%              a=0.001;
%              colr = ['-or'; '-ob'; '-om';'-og';'-oc'];
%              for ISEG=1:NSEG
%                  XYZSEG(:,NOD1)=XYZVERT(:,SEG_CELL(NOD1,ISEG));
%                  XYZSEG(:,NOD2)=XYZVERT(:,SEG_CELL(NOD2,ISEG));
%                  IG         =max(1,SEG_CELL(6,ISEG));
%                  plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),colr(IG,1:3),'LineWidth',2)
%                  text(0.5*(sum(XYZSEG(IAXIS,:)))+a,0.5*(sum(XYZSEG(JAXIS,:)))+a,...
%                       0.5*(sum(XYZSEG(KAXIS,:)))+a,num2str(SEG_CELL(4,ISEG)),'FontSize',10)
%              end
%              for IVERT=1:NVERT
%                  text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                       XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
%
%              end
%              xlabel('X')
%              ylabel('Y')
%              zlabel('Z')
%              %axis equal
%              box on
%              view([45 45])
%              pause
%          end


         if(NSEG < 3 ); continue; end

         % Ear clipping algorithm by TRIANGLE and BODY:
         % 1. Define closed 3D polyline:
         [IFLG,NPOLY,ILO_POLY,NSG_POLY,NSEG,SEG_CELL,SEG_POS]=GET_CLOSED_POLYLINES(NSEG,SEG_CELL,SEG_POS);


         % Figure:
         if (plotflg || (IFLG ~=0))
         scrsz = get(groot,'ScreenSize');
         figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%          subplot(1,3,1)
         hold on
         colr = ['-or'; '-ob'; '-om'];
%          for ISEG=1:NSEG
%              XYZSEG(:,1)=XYZVERT(:,SEG_CELL(1,ISEG));
%              XYZSEG(:,2)=XYZVERT(:,SEG_CELL(2,ISEG));
%              IG         =SEG_CELL(6,ISEG);
%              plot3(XYZSEG(1,:),XYZSEG(2,:),XYZSEG(3,:),colr(IG,1:3),'LineWidth',2)
%
%              if(BNDINT_FLAG)
%              P(1,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J-1) MESHES(NM).Z(K-1)];
%              P(2,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J-1) MESHES(NM).Z(K-1)];
%              P(3,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J  ) MESHES(NM).Z(K-1)];
%              P(4,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J  ) MESHES(NM).Z(K-1)];
%              P(5,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J-1) MESHES(NM).Z(K  )];
%              P(6,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J-1) MESHES(NM).Z(K  )];
%              P(7,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J  ) MESHES(NM).Z(K  )];
%              P(8,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J  ) MESHES(NM).Z(K  )];
%
%              FC(1,:) = [ 1 4 8 5 ];
%              FC(2,:) = [ 2 3 7 6 ];
%              FC(3,:) = [ 1 5 6 2 ];
%              FC(4,:) = [ 4 8 7 3 ];
%              FC(5,:) = [ 1 2 3 4 ];
%              FC(6,:) = [ 5 6 7 8 ];
%
%              for IFC=1:6
%                  [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
%                  set(hp,'FaceAlpha',0.1)
%              end
%              end
%
%          end
%          xlabel('X')
%          ylabel('Y')
%          zlabel('Z')
%          axis equal
%          box on
%          view([45 45])
         end

         if(IFLG ~=0)
             disp(['IFLG ~=0, couldn''t close polyline, ' num2str(BNDINT_FLAG) ': ' ...
                    num2str(I) ', ' num2str(J) ', ' num2str(K)])
             NSEG
             SEG_CELL
             XYZVERT
         end
         if (plotflg || (IFLG ~=0))
         a=0.02/3*(DXCELL(I)+DYCELL(J)+DZCELL(K));
         subplot(1,2,1)
         hold on
         [ierr]=plot_geometries(N_GEOMETRY,GEOM,0);
         colr = ['-or'; '-ob'; '-om'];
         for ISEG=1:NSEG
             XYZSEG(:,1)=XYZVERT(:,SEG_CELL(1,ISEG));
             XYZSEG(:,2)=XYZVERT(:,SEG_CELL(2,ISEG));
             IG         =max(1,SEG_CELL(6,ISEG));
             plot3(XYZSEG(1,:),XYZSEG(2,:),XYZSEG(3,:),colr(IG,1:3),'LineWidth',2)

             text(0.5*(sum(XYZSEG(IAXIS,:)))+a,0.5*(sum(XYZSEG(JAXIS,:)))+a,...
                  0.5*(sum(XYZSEG(KAXIS,:)))+a,num2str(SEG_CELL(4,ISEG)),'FontSize',10,'Color','b')

             if(BNDINT_FLAG)
             P(1,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J-1) MESHES(NM).Z(K-1)];
             P(2,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J-1) MESHES(NM).Z(K-1)];
             P(3,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J  ) MESHES(NM).Z(K-1)];
             P(4,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J  ) MESHES(NM).Z(K-1)];
             P(5,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J-1) MESHES(NM).Z(K  )];
             P(6,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J-1) MESHES(NM).Z(K  )];
             P(7,IAXIS:KAXIS) = [MESHES(NM).X(I  ) MESHES(NM).Y(J  ) MESHES(NM).Z(K  )];
             P(8,IAXIS:KAXIS) = [MESHES(NM).X(I-1) MESHES(NM).Y(J  ) MESHES(NM).Z(K  )];

             FC(1,:) = [ 1 4 8 5 ];
             FC(2,:) = [ 2 3 7 6 ];
             FC(3,:) = [ 1 5 6 2 ];
             FC(4,:) = [ 4 8 7 3 ];
             FC(5,:) = [ 1 2 3 4 ];
             FC(6,:) = [ 5 6 7 8 ];

             for IFC=1:6
                 [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
                 set(hp,'FaceAlpha',0.1)
             end
             end

         end
         for IVERT=1:NVERT
             text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
                  XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10,'Color','r')

         end
         title(['I,J,K=' num2str(I) ', ' num2str(J) ', ' num2str(K)])
         xlabel('X')
         ylabel('Y')
         zlabel('Z')
         axis equal
         box on
         view([45 45])
         subplot(1,2,2)
         hold on
         for FAXIS=IAXIS:KAXIS
            for ILH=LOW_IND:HIGH_IND
               switch(FAXIS)
                   case(IAXIS)
                       IADD=ILH-2;  JADD=    0;  KADD=    0;
                   case(JAXIS)
                       IADD=    0;  JADD=ILH-2;  KADD=    0;
                   case(KAXIS)
                       IADD=    0;  JADD=    0;  KADD=ILH-2;
               end

               % By segment: Add Vertices/Segments to local arrays:
               ICF = MESHES(NM).FCVAR(I+IADD,J+JADD,K+KADD,IBM_IDCF,FAXIS); %CEIB_XYZ(ILH,FAXIS);
               if ( ICF > 0 )
                   NFACE = MESHES(NM).CUT_FACE(ICF).NFACE;
                   XYZV  = MESHES(NM).CUT_FACE(ICF).XYZVERT;
                   for JCF=1:NFACE
                       NELEM  = MESHES(NM).CUT_FACE(ICF).CFELEM(1,JCF);
                       CFELEM = MESHES(NM).CUT_FACE(ICF).CFELEM(2:NELEM+1,JCF);
                       [hp]=patch(XYZV(IAXIS,CFELEM),XYZV(JAXIS,CFELEM),XYZV(KAXIS,CFELEM),'r');
                       set(hp,'FaceAlpha',0.3)
                   end
               end
            end
         end
         xlabel('X')
         ylabel('Y')
         zlabel('Z')
         axis equal
         box on
         end

         % 2. Define triangles by Body and triangle, all triangles generated
         %    point outside of solid region:
         [NFACE, CFELEM, BOD_TRI] = EAR_CLIP(NSEG,NVERT,SEG_CELL,XYZVERT,I,J,K,NPOLY,ILO_POLY,NSG_POLY);


         if(plotflg || (IFLG ~=0))
             % Plot faces:
             for ICF=1:NFACE
                 XYZFC(:,NOD1)=XYZVERT(:,CFELEM(2,ICF));
                 XYZFC(:,NOD2)=XYZVERT(:,CFELEM(3,ICF));
                 XYZFC(:,NOD3)=XYZVERT(:,CFELEM(4,ICF));
                 IG         =BOD_TRI(1,ICF);
                 [hp]=patch(XYZFC(IAXIS,:),XYZFC(JAXIS,:),XYZFC(KAXIS,:),colr(IG,1:3));
                 set(hp,'FaceAlpha',0.1)
                 text(1/3*(sum(XYZFC(IAXIS,:)))+a,1/3*(sum(XYZFC(JAXIS,:)))+a,...
                      1/3*(sum(XYZFC(KAXIS,:)))+a,num2str(BOD_TRI(2,ICF)),'FontSize',10)
             end
             view([45 45])
             pause
             close
         end

         % This is a cut-face, allocate space:
         NCUTFACE = MESHES(NM).N_CUTFACE_MESH + MESHES(NM).N_GCCUTFACE_MESH + 1;
         if (BNDINT_FLAG)
            MESHES(NM).N_CUTFACE_MESH = NCUTFACE;
         else
            MESHES(NM).N_GCCUTFACE_MESH = MESHES(NM).N_GCCUTFACE_MESH + 1;
         end
         MESHES(NM).CCVAR(I,J,K,IBM_IDCF) = NCUTFACE;
         MESHES(NM).CUT_FACE(NCUTFACE).NVERT  = NVERT;
         MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = 0;
         MESHES(NM).CUT_FACE(NCUTFACE).IJK(1:MAX_DIM+1) = [ I, J, K, 0 ]; % No axis = 0
         MESHES(NM).CUT_FACE(NCUTFACE).STATUS = IBM_INBOUNDARY;
         MESHES(NM).CUT_FACE(NCUTFACE).XYZVERT(IAXIS:KAXIS,1:NVERT) = XYZVERT(IAXIS:KAXIS,1:NVERT);

         % Assign surf-index: Depending on GEOMETRY:
         NCF = 0;
         for ICF=1:NFACE
            IBOD = BOD_TRI(1,ICF); ITRI = BOD_TRI(2,ICF);

            % Area properties for special cfaces:
            % Computed from the cross product:
            D23 = XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD3,ICF)) - ...
                  XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD2,ICF));
            D12 = XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD2,ICF)) - ...
                  XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD1,ICF));
            NORMTRI(IAXIS:KAXIS) = cross(D12',D23');
            NNORM                = norm(NORMTRI);
            if (NNORM < 2.*GEOMEPS^2.); continue; end
            NORMTRI(IAXIS:KAXIS) = NORMTRI(IAXIS:KAXIS) / NNORM;

            % First test if INB face is on Cartesian face and pointing
            % outside of Cartesian cell. If so drop:
            XYZ(IAXIS:KAXIS) = XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD1,ICF))';
            % IAXIS:
            if ( (abs(NORMTRI(IAXIS)+1.) < GEOMEPS) && ...
                 (abs(XFACE(I-FCELL  )-XYZ(IAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(IAXIS)-1.) < GEOMEPS) && ...
                 (abs(XFACE(I-FCELL+1)-XYZ(IAXIS)) < GEOMEPS) ); continue; end % High Face
            % JAXIS:
            if ( (abs(NORMTRI(JAXIS)+1.) < GEOMEPS) && ...
                 (abs(YFACE(J-FCELL  )-XYZ(JAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(JAXIS)-1.) < GEOMEPS) && ...
                 (abs(YFACE(J-FCELL+1)-XYZ(JAXIS)) < GEOMEPS) ); continue; end % High Face
            % KAXIS:
            if ( (abs(NORMTRI(KAXIS)+1.) < GEOMEPS) && ...
                 (abs(ZFACE(K-FCELL  )-XYZ(KAXIS)) < GEOMEPS) ); continue; end % Low Face
            if ( (abs(NORMTRI(KAXIS)-1.) < GEOMEPS) && ...
                 (abs(ZFACE(K-FCELL+1)-XYZ(KAXIS)) < GEOMEPS) ); continue; end % High Face

            % Area:
            AREA  = 0.5*NNORM;

            % Face Vertices average location:
            ACEN(IAXIS:KAXIS) = 1./3.*(XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD1,ICF))' + ...
                                       XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD2,ICF))' + ...
                                       XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD3,ICF))');

            % dot(i,nc) int(x)dA
            INXAREA = NORMTRI(IAXIS)*ACEN(IAXIS)*AREA; % Single Gauss pt integration.

            XC1(IAXIS:KAXIS) = 0.5*(XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD2,ICF)) + ...
                                    XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD3,ICF)))'; % X23
            XC2(IAXIS:KAXIS) = 0.5*(XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD1,ICF)) + ...
                                    XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD3,ICF)))';
            X12(IAXIS:KAXIS) = 0.5*(XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD1,ICF)) + ...
                                    XYZVERT(IAXIS:KAXIS,CFELEM(1+NOD2,ICF)))';
            % dot(i,nc) int(x^2)dA, dot(j,nc) int(y^2)dA, dot(k,nc) int(z^2)dA
            SQAREA(IAXIS:KAXIS) = 0.;
            for IX=IAXIS:KAXIS
                INT2 = (XC1(IX)^2. + XC2(IX)^2. + X12(IX)^2.) / 3.;
                SQAREA(IX) = SQAREA(IX) + NORMTRI(IX)*INT2*AREA;  % Midpoint rule.
            end

            NCF = NCF + 1;

            MESHES(NM).CUT_FACE(NCUTFACE).AREA(NCF) = AREA;
            MESHES(NM).CUT_FACE(NCUTFACE).XYZCEN(IAXIS:KAXIS,NCF) = ACEN(IAXIS:KAXIS)';

            % Fields for cut-cell volume/centroid computation:
            % dot(i,nc)*int(x)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).INXAREA(NCF)   = INXAREA;
            % dot(i,nc)*int(x^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).INXSQAREA(NCF) = SQAREA(IAXIS);
            % dot(j,nc)*int(y^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).JNYSQAREA(NCF) = SQAREA(JAXIS);
            % dot(k,nc)*int(z^2)dA:
            MESHES(NM).CUT_FACE(NCUTFACE).KNZSQAREA(NCF) = SQAREA(KAXIS);

            % Define Body-triangle reference:
            MESHES(NM).CUT_FACE(NCUTFACE).BODTRI(1:2,NCF)= [IBOD ITRI]';

            MESHES(NM).CUT_FACE(NCUTFACE).SURF_INDEX(NCF) = GEOM(IBOD).SURFS(ITRI);

            % All faces connectivities:
            MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1:1+NOD3,NCF) = CFELEM(1:1+NOD3,ICF);

         end

         MESHES(NM).CUT_FACE(NCUTFACE).NFACE = NCF;


      end % I
   end % J
end % K


ierr=0;

return
