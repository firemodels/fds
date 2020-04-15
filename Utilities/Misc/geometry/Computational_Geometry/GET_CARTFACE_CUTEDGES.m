function [ierr]=GET_CARTFACE_CUTEDGES(X1AXIS,X2AXIS,X3AXIS,             ...
                                      XIAXIS,XJAXIS,XKAXIS,NM,          ...
                                      X2LO,X2HI,X3LO,X3HI,X2LO_CELL,X2HI_CELL, ...
                                      X3LO_CELL,X3HI_CELL,INDX1,X1PLN)

global IAXIS JAXIS KAXIS NOD1 NOD2 IBM_MAX_WSTRIANG_SEG LOW_IND HIGH_IND FCELL 
global IBM_INBOUNDCF IBM_GG IBM_IDCE                             
global MESHES BODINT_PLANE
global GEOMEPS MAX_DIM
global X2NOC X3NOC
global X1FACE X2FACE DX2FACE X3FACE DX3FACE
global FACERT CELLRT

ierr=1;

% Segment by segment define the INBOUNDARY MESHES(NM).IBM_CUT_EDGES between crossings
% and individualize the Cartesian face they belong to.    
for ISEG=1:BODINT_PLANE.NSEGS

   SEG(NOD1:NOD2)    = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD1));
   XYZ2(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD2));

   if (max(XYZ1(X2AXIS),XYZ2(X2AXIS)) < X2FACE(X2LO)-GEOMEPS); continue; end
   if (min(XYZ1(X2AXIS),XYZ2(X2AXIS)) > X2FACE(X2HI)+GEOMEPS); continue; end
   if (max(XYZ1(X3AXIS),XYZ2(X3AXIS)) < X3FACE(X3LO)-GEOMEPS); continue; end
   if (min(XYZ1(X3AXIS),XYZ2(X3AXIS)) > X3FACE(X3HI)+GEOMEPS); continue; end
   
   NBCROSS = BODINT_PLANE.NBCROSS(ISEG); % Cross points include Node1, Node2

   
   % x2_x3 of segment point 1:
   X2_1 = XYZ1(X2AXIS); X3_1 = XYZ1(X3AXIS);
   % x2_x3 of segment point 2:
   X2_2 = XYZ2(X2AXIS); X3_2 = XYZ2(X3AXIS);
   
   % Normal out:
   SLEN = sqrt( (X2_2-X2_1)^2. + (X3_2-X3_1)^2. );
   STANI(IAXIS:JAXIS) = [ (X2_2-X2_1), (X3_2-X3_1) ]*SLEN^(-1.);
   SNORI(IAXIS:JAXIS) = [ STANI(JAXIS), -STANI(IAXIS) ];

   INDSEG(1:IBM_MAX_WSTRIANG_SEG+2) = BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2, ISEG);
   NTRISEG = INDSEG(1);

   ADD2FACES = false;
   % Type to be assigned to cut edges:
   CETYPE = 2*(BODINT_PLANE.SEGTYPE(LOW_IND,ISEG)+1) - BODINT_PLANE.SEGTYPE(HIGH_IND,ISEG);
   if ( CETYPE == IBM_GG ); ADD2FACES = true; end

   INRAY  = false;

   % Different cases:
   % First check if segment geomepsilon aligned with x2:
   if (BODINT_PLANE.X2ALIGNED(ISEG))

      % Test if node1 of segment is in geomepsilon vicinity of an x2 ray
      for KK=X3LO:X3HI
         % x3 location of ray along x2, on the x2-x3 plane:
         X3RAY = X3FACE(KK);
         if ( abs(X3RAY-X3_1) < GEOMEPS )
            INRAY = true;
            break
         end
      end

      if (INRAY) % Segment in x2 ray defined by x3 face index kk.

         % 1. INB cut-edges on top of an x2 gridline, assign to cut-face
         %    defined by normal out.
         KK2VEC(LOW_IND:HIGH_IND) = 0;
         if (ADD2FACES)
             NPFACE   = 2;
             KK2VEC(LOW_IND) = KK + FCELL;
             KK2VEC(HIGH_IND)= KK + FCELL - 1;
         else
             NPFACE = 1;
             if ( SNORI(JAXIS) > 0. ) % add 1 to index kk+FCELL-1 (i.e. lower face index)
                 KK2VEC(LOW_IND) = KK + FCELL;
             else
                 KK2VEC(LOW_IND)= KK + FCELL - 1;
             end
         end

         for IPFACE=1:NPFACE

            KK2 = KK2VEC(IPFACE);

            % Figure out which cut faces the inboundary cut-edges of
            % this segment belong to:
            % We have nbcross-1 INBOUNDARY CUT_EDGEs to generate.
            for IEDGE=1:NBCROSS-1

               % Location along Segment:
               SVAR1 = BODINT_PLANE.SVAR(IEDGE  ,ISEG);
               SVAR2 = BODINT_PLANE.SVAR(IEDGE+1,ISEG);
               % Location of midpoint of cut-edge:
               SVAR12 = 0.5 * (SVAR1+SVAR2);
              
               % Define Cartesian segment this cut-edge belongs:
               XPOS   = X2_1 + SVAR12*STANI(IAXIS);
               if (X2NOC==0)
                  JJ2 = floor((XPOS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
                  % Discard cut-edges on faces laying on x2hi.
                  if ((JJ2 < X2LO_CELL) || (JJ2 > X2HI_CELL)); continue; end
               else
                  FOUND_SEG=false;
                  for JJ2=X2LO_CELL:X2HI_CELL
                     % Check if XPOS is within this segment JJ2:
                     if ((XPOS-X2FACE(JJ2-1)) >= 0. && (X2FACE(JJ2)-XPOS) > 0.)
                        FOUND_SEG=true;
                        break
                     end
                  end
                  if (~FOUND_SEG); continue; end
               end

               if ((KK2 < X3LO_CELL) || (KK2 > X3HI_CELL)); continue; end
               
               % HERE IF NEEDED TEST IF SEG IS INSIDE OR OUTSIDE.
               % If segment is inside the solid region mark cells surrounding face
               % to be treated in special manner (only if they happen to be type CUTCFE), 
               % then drop segment.
               XY(IAXIS:JAXIS) = [X2_1 X3_1] + SVAR12*STANI(IAXIS:JAXIS);
               [IS_SOLID] = GET_IS_SOLID_PT(BODINT_PLANE,X1AXIS,X2AXIS,X3AXIS,XY,SNORI,X1PLN);
               if (IS_SOLID); continue; end
               
               % Face indexes:
               INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ2, KK2 ]; % Local x1,x2,x3
               INDIF=INDXI(XIAXIS);
               INDJF=INDXI(XJAXIS);
               INDKF=INDXI(XKAXIS);           
                              
               % Now the face is, FCVAR (x1axis):
               if (MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS) > 0) % There is already
                                                                                 % an entry in CUT_EDGE.
                  CEI = MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS);
               else % We need a new entry in CUT_EDGE
                  CEI      = MESHES(NM).N_CUTEDGE_MESH + 1;
                  MESHES(NM).N_CUTEDGE_MESH                       = CEI;
                  MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS)= CEI;
                  %CALL CUT_EDGE_ARRAY_REALLOC(NM,CEI)
                  MESHES(NM).CUT_EDGE(CEI).NVERT   = 0;
                  %CALL NEW_EDGE_ALLOC(NM,CEI,IBM_ALLOC_DVERT,IBM_ALLOC_DELEM)
                  MESHES(NM).CUT_EDGE(CEI).NEDGE   = 0;
                  MESHES(NM).CUT_EDGE(CEI).NEDGE1  = 0;
                  MESHES(NM).CUT_EDGE(CEI).IJK(1:MAX_DIM+2) = [ INDIF, INDJF, INDKF, X1AXIS, CETYPE ];
                  MESHES(NM).CUT_EDGE(CEI).STATUS  = IBM_INBOUNDCF;
               end

               % Add vertices, non repeated vertex entries at this point.
               NVERT = MESHES(NM).CUT_EDGE(CEI).NVERT;
               % Define vertices for this segment:
               %                                           xv1                      yv1                     zv1
               XYZV1LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR1, X3_1+STANI(JAXIS)*SVAR1 ];
               XYZV1(IAXIS) = XYZV1LC(XIAXIS);
               XYZV1(JAXIS) = XYZV1LC(XJAXIS);
               XYZV1(KAXIS) = XYZV1LC(XKAXIS);
               [NVERT,INOD1]=INSERT_FACE_VERT(XYZV1,NM,CEI,NVERT);
               %                                           xv2                      yv2                     zv2
               XYZV2LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR2, X3_1+STANI(JAXIS)*SVAR2 ];
               XYZV2(IAXIS) = XYZV2LC(XIAXIS);
               XYZV2(JAXIS) = XYZV2LC(XJAXIS);
               XYZV2(KAXIS) = XYZV2LC(XKAXIS);
               [NVERT,INOD2]=INSERT_FACE_VERT(XYZV2,NM,CEI,NVERT);

               NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
               %CALL REALLOCATE_EDGE_ELEM(NM,CEI,NEDGE+1)
               if ( NPFACE == 1 )
                  MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
               else
                  DIRAXIS = X2AXIS;
                  CONDAX  = (XYZV2(DIRAXIS)-XYZV1(DIRAXIS)) > 0;
                  if ( KK2 == KK+FCELL-1 )
                     if (CONDAX)
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
                     else
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD2, INOD1 ];
                     end
                  else
                     if (CONDAX)
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD2, INOD1 ];
                     else
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
                     end
                  end
               end
               MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,NEDGE+1) = ...
                               BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,ISEG);
               MESHES(NM).CUT_EDGE(CEI).INDSEG(IBM_MAX_WSTRIANG_SEG+3,NEDGE+1) = ...
                                     -sum(BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG))/2;
               MESHES(NM).CUT_EDGE(CEI).NVERT = NVERT;
               MESHES(NM).CUT_EDGE(CEI).NEDGE = NEDGE+1;
               MESHES(NM).CUT_EDGE(CEI).NEDGE1= MESHES(NM).CUT_EDGE(CEI).NEDGE;
               
               % Test for Repeated edge -> If so note FACERT
               for IDG=1:NEDGE
                   if( ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD1   && ...
                         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD2 ) || ...
                       ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD2   && ...
                         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD1 ) )
                       FACERT(JJ2,KK2) = 1;
                       break
                   end
               end
      
            end
         end
         continue % Skips rest of iseg loop, for this ISEG.
      end

   % Second check if segment geomepsilon aligned with x3:
   elseif (BODINT_PLANE.X3ALIGNED(ISEG))

      % Test if node1 of segment is in geomepsilon vicinity of an x3 ray
      for JJ=X2LO:X2HI
         % x2 location of ray along x3, on the x2-x3 plane:
         X2RAY = X2FACE(JJ);
         if ( abs(X2RAY-X2_1) < GEOMEPS )
            INRAY = true;
            break
         end
      end

      if (INRAY) % Segment in x3 ray defined by x2 face index JJ

         % 1. INB cut-edges on top of an x3 gridline, assign to cut-face
         %    defined by normal out.
         JJ2VEC(LOW_IND:HIGH_IND) = 0;
         if (ADD2FACES)
            NPFACE = 2;
            JJ2VEC(LOW_IND)  = JJ + FCELL;
            JJ2VEC(HIGH_IND) = JJ + FCELL - 1;
         else
            NPFACE = 1;
            if ( SNORI(IAXIS) > 0. ) % add 1 to index jj+FCELL-1 (i.e. lower face index)
               JJ2VEC(LOW_IND) = JJ + FCELL;
            else
               JJ2VEC(LOW_IND) = JJ + FCELL - 1;
            end
         end

         for IPFACE=1:NPFACE

            JJ2 = JJ2VEC(IPFACE);

            % Figure out which cut faces the inboundary cut-edges of
            % this segment belong to:
            % We have NBCROSS-1 INBOUNDARY CUT_EDGEs to generate.
            for IEDGE=1:NBCROSS-1

               % Location along Segment:
               SVAR1 = BODINT_PLANE.SVAR(IEDGE  ,ISEG);
               SVAR2 = BODINT_PLANE.SVAR(IEDGE+1,ISEG);
               % Location of midpoint of cut-edge:
               SVAR12 = 0.5 * (SVAR1+SVAR2);               
               
               % Define Cartesian segment this cut-edge belongs:
               XPOS = X3_1 + SVAR12*STANI(JAXIS);
               if (X3NOC==0)
                  KK2 = floor((XPOS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
                  % Discard cut-edges on faces laying on x3hi.
                  if ((KK2 < X3LO_CELL) || (KK2 > X3HI_CELL)); continue; end
               else
                  FOUND_SEG=false;
                  for KK2=X3LO_CELL:X3HI_CELL
                     % Check if XPOS is within this segment KK2:
                     if ((XPOS-X3FACE(KK2-1)) >= 0. && (X3FACE(KK2)-XPOS) > 0.)
                        FOUND_SEG=true;
                        break
                     end
                  end
                  if (~FOUND_SEG); continue; end
               end

               if ((JJ2 < X2LO_CELL) || (JJ2 > X2HI_CELL)); continue; end
               
               % HERE IF NEEDED TEST IF SEG IS INSIDE OR OUTSIDE.
               % If segment is inside the solid region mark cells surrounding face
               % to be treated in special manner (only if they happen to be type CUTCFE), 
               % then drop segment.
               XY(IAXIS:JAXIS) = [X2_1 X3_1] + SVAR12*STANI(IAXIS:JAXIS);
               [IS_SOLID] = GET_IS_SOLID_PT(BODINT_PLANE,X1AXIS,X2AXIS,X3AXIS,XY,SNORI,X1PLN);
               if (IS_SOLID); continue; end
               
               % Face indexes:
               INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ2, KK2 ]; % Local x1,x2,x3
               INDIF=INDXI(XIAXIS);
               INDJF=INDXI(XJAXIS);
               INDKF=INDXI(XKAXIS);
                              
               % Now the face is, FCVAR (x1axis):
               if (MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS) > 0) % There is already
                                                                            % an entry in CUT_EDGE.
                  CEI = MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS);
               else % We need a new entry in CUT_EDGE
                  CEI      = MESHES(NM).N_CUTEDGE_MESH + 1;
                  MESHES(NM).N_CUTEDGE_MESH                       = CEI;
                  MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS)= CEI;
                  %CALL CUT_EDGE_ARRAY_REALLOC(NM,CEI)
                  MESHES(NM).CUT_EDGE(CEI).NVERT   = 0;
                  %CALL NEW_EDGE_ALLOC(NM,CEI,IBM_ALLOC_DVERT,IBM_ALLOC_DELEM)
                  MESHES(NM).CUT_EDGE(CEI).NEDGE   = 0;
                  MESHES(NM).CUT_EDGE(CEI).NEDGE1  = 0;
                  MESHES(NM).CUT_EDGE(CEI).IJK(1:MAX_DIM+2) = [ INDIF, INDJF, INDKF, X1AXIS, CETYPE ];
                  MESHES(NM).CUT_EDGE(CEI).STATUS  = IBM_INBOUNDCF;
               end

               % Add vertices, non repeated vertex entries at this point.
               NVERT = MESHES(NM).CUT_EDGE(CEI).NVERT;
               % Define vertices for this segment:
               %                                           xv1                      yv1                     zv1
               XYZV1LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR1, X3_1+STANI(JAXIS)*SVAR1 ];
               XYZV1(IAXIS) = XYZV1LC(XIAXIS);
               XYZV1(JAXIS) = XYZV1LC(XJAXIS);
               XYZV1(KAXIS) = XYZV1LC(XKAXIS);
               [NVERT,INOD1]=INSERT_FACE_VERT(XYZV1,NM,CEI,NVERT);
               %                                           xv2                      yv2                     zv2
               XYZV2LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR2, X3_1+STANI(JAXIS)*SVAR2 ];
               XYZV2(IAXIS) = XYZV2LC(XIAXIS);
               XYZV2(JAXIS) = XYZV2LC(XJAXIS);
               XYZV2(KAXIS) = XYZV2LC(XKAXIS);
               [NVERT,INOD2]=INSERT_FACE_VERT(XYZV2,NM,CEI,NVERT);

               NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
               %CALL REALLOCATE_EDGE_ELEM(NM,CEI,NEDGE+1)
               if ( NPFACE == 1 )
                  MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
               else
                  DIRAXIS = X3AXIS;
                  CONDAX  = (XYZV2(DIRAXIS)-XYZV1(DIRAXIS)) > 0;
                  if ( JJ2 == JJ+FCELL-1 )
                     if (CONDAX)
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD2, INOD1 ];
                     else
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
                     end
                  else
                     if (CONDAX)
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
                     else
                        MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD2, INOD1 ];
                     end
                  end
               end
               MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,NEDGE+1) = ...
                               BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,ISEG);
               MESHES(NM).CUT_EDGE(CEI).INDSEG(IBM_MAX_WSTRIANG_SEG+3,NEDGE+1) = ...
                                     -sum(BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG))/2;
               MESHES(NM).CUT_EDGE(CEI).NVERT = NVERT;
               MESHES(NM).CUT_EDGE(CEI).NEDGE = NEDGE+1;
               MESHES(NM).CUT_EDGE(CEI).NEDGE1= MESHES(NM).CUT_EDGE(CEI).NEDGE;

               % Test for Repeated edge -> If so note FACERT
               for IDG=1:NEDGE
                   if( ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD1   && ...
                         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD2 ) || ...
                       ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD2   && ...
                         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD1 ) )
                       FACERT(JJ2,KK2) = 1;
                       break
                   end
               end
               
            end
         end
         continue % Skips rest of iseg loop, for this ISEG.
      end

   end

   % 3. Regular case: INB cut-edge with centroid inside a
   %    Cartesian face, assign to corresponding FCVAR IBM_IDCE variable.
   % This is the most common case, INBOUNDARY edges defined inside x1 faces.
   % We have NBCROSS-1 INBOUNDARY CUT_EDGEs to generate.
   for IEDGE=1:NBCROSS-1

      % Location along Segment:
      SVAR1 = BODINT_PLANE.SVAR(IEDGE  ,ISEG);
      SVAR2 = BODINT_PLANE.SVAR(IEDGE+1,ISEG);
      % Location of midpoint of cut-edge:
      SVAR12 = 0.5 * (SVAR1+SVAR2);
      
      % Define Cartesian face this cut-edge belongs:
      XPOS = X2_1 + SVAR12*STANI(IAXIS);
      if (X2NOC==0)
         JJ2  = floor((XPOS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
         if ((JJ2 < X2LO_CELL) || (JJ2 > X2HI_CELL)); continue; end
      else
         FOUND_SEG=false;
         for JJ2=X2LO_CELL:X2HI_CELL
            % Check if XPOS is within this segment JJ2:
            if ((XPOS-X2FACE(JJ2-1)) >= 0. && (X2FACE(JJ2)-XPOS) > 0.)
               FOUND_SEG=true;
               break
            end
         end
         if (~FOUND_SEG); continue; end
      end
      XPOS = X3_1 + SVAR12*STANI(JAXIS);
      if (X3NOC==0)
         KK2  = floor((XPOS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
         if ((KK2 < X3LO_CELL) || (KK2 > X3HI_CELL)); continue; end
      else
         FOUND_SEG=false;
         for KK2=X3LO_CELL:X3HI_CELL
            % Check if XPOS is within this segment KK2:
            if ((XPOS-X3FACE(KK2-1)) >= 0. && (X3FACE(KK2)-XPOS) > 0.)
               FOUND_SEG=true;
               break
            end
         end
         if (~FOUND_SEG); continue; end
      end
   
      % HERE IF NEEDED TEST IF SEG IS INSIDE OR OUTSIDE.
      % If segment is inside the solid region mark cells surrounding face
      % to be treated in special manner (only if they happen to be type CUTCFE),
      % then drop segment.
      XY(IAXIS:JAXIS) = [X2_1 X3_1] + SVAR12*STANI(IAXIS:JAXIS);
      [IS_SOLID] = GET_IS_SOLID_PT(BODINT_PLANE,X1AXIS,X2AXIS,X3AXIS,XY,SNORI,X1PLN);
      if (IS_SOLID); continue; end
      
      % Face indexes:
      INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ2, KK2 ]; % Local x1,x2,x3
      INDIF=INDXI(XIAXIS);
      INDJF=INDXI(XJAXIS);
      INDKF=INDXI(XKAXIS);
      
      % Now the face is, FCVAR (x1axis):
      if (MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS) > 0) % There is already
                                                                   % an entry in CUT_EDGE.
         CEI = MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS);
      else % We need a new entry in CUT_EDGE
         CEI      = MESHES(NM).N_CUTEDGE_MESH + 1;
         MESHES(NM).N_CUTEDGE_MESH                       = CEI;
         MESHES(NM).FCVAR(INDIF,INDJF,INDKF,IBM_IDCE,X1AXIS)= CEI;
         %CALL CUT_EDGE_ARRAY_REALLOC(NM,CEI)
         MESHES(NM).CUT_EDGE(CEI).NVERT   = 0;
         %CALL NEW_EDGE_ALLOC(NM,CEI,IBM_ALLOC_DVERT,IBM_ALLOC_DELEM)
         MESHES(NM).CUT_EDGE(CEI).NEDGE   = 0;
         MESHES(NM).CUT_EDGE(CEI).NEDGE1  = 0;
         MESHES(NM).CUT_EDGE(CEI).IJK(1:MAX_DIM+2) = [ INDIF, INDJF, INDKF, X1AXIS, CETYPE ];
         MESHES(NM).CUT_EDGE(CEI).STATUS  = IBM_INBOUNDCF;
      end

      % Add vertices, non repeated vertex entries at this point.
      NVERT = MESHES(NM).CUT_EDGE(CEI).NVERT;

      % Define vertices for this segment:
      %                                           xv1                      yv1                     zv1
      XYZV1LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR1, X3_1+STANI(JAXIS)*SVAR1 ];
      XYZV1(IAXIS) = XYZV1LC(XIAXIS);
      XYZV1(JAXIS) = XYZV1LC(XJAXIS);
      XYZV1(KAXIS) = XYZV1LC(XKAXIS);
      [NVERT,INOD1]=INSERT_FACE_VERT(XYZV1,NM,CEI,NVERT);
      %                                           xv2                      yv2                     zv2
      XYZV2LC(IAXIS:KAXIS)= [ X1FACE(INDX1(X1AXIS)), X2_1+STANI(IAXIS)*SVAR2, X3_1+STANI(JAXIS)*SVAR2 ];
      XYZV2(IAXIS) = XYZV2LC(XIAXIS);
      XYZV2(JAXIS) = XYZV2LC(XJAXIS);
      XYZV2(KAXIS) = XYZV2LC(XKAXIS);
      [NVERT,INOD2]=INSERT_FACE_VERT(XYZV2,NM,CEI,NVERT);

      NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
      %CALL REALLOCATE_EDGE_ELEM(NM,CEI,NEDGE+1)
      MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE+1) = [ INOD1, INOD2 ];
      MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,NEDGE+1) = ...
                      BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,ISEG);
      MESHES(NM).CUT_EDGE(CEI).INDSEG(IBM_MAX_WSTRIANG_SEG+3,NEDGE+1) = ...
                            -sum(BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG))/2;
      MESHES(NM).CUT_EDGE(CEI).NVERT = NVERT;
      MESHES(NM).CUT_EDGE(CEI).NEDGE = NEDGE+1;
      MESHES(NM).CUT_EDGE(CEI).NEDGE1= MESHES(NM).CUT_EDGE(CEI).NEDGE;
      
      % Test for Repeated edge -> If so note FACERT
      for IDG=1:NEDGE
          if( ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD1   && ...
                MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD2 ) || ...
              ( MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1,IDG) == INOD2   && ...
                MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD2,IDG) == INOD1 ) )
              FACERT(JJ2,KK2) = 1;
              break
          end
      end
   end

end

% Note cells in CELLRT due to FCERT intersections in GET_BODINT_PLANE:
for KK2=X3LO_CELL:X3HI_CELL
    for JJ2=X2LO_CELL:X2HI_CELL 
       if(~FACERT(JJ2,KK2)); continue; end
       % Low cell indexes:
       INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ2, KK2 ]; % Local x1,x2,x3
       INDIF=INDXI(XIAXIS);
       INDJF=INDXI(XJAXIS);
       INDKF=INDXI(XKAXIS);
       
       CELLRT(INDIF,INDJF,INDKF) = true;
       
       % High cell indexes:
       INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS)+1, JJ2, KK2 ]; % Local x1,x2,x3
       INDIF=INDXI(XIAXIS);
       INDJF=INDXI(XJAXIS);
       INDKF=INDXI(XKAXIS);
       
       CELLRT(INDIF,INDJF,INDKF) = true;           

    end
end

ierr=0;

return