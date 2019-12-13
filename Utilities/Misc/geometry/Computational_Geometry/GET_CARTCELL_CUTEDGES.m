function [ierr]=GET_CARTCELL_CUTEDGES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND)

global N_GEOMETRY GEOM IBM_DELTA_NBCROSS NOD1 NOD2 IAXIS JAXIS KAXIS
global MAX_DIM CCGUARD GEOMEPS TRANS MESHES
global XFACE DXFACE YFACE DYFACE ZFACE DZFACE
global ILO_FACE IHI_FACE JLO_FACE JHI_FACE KLO_FACE KHI_FACE
global ILO_CELL IHI_CELL JLO_CELL JHI_CELL KLO_CELL KHI_CELL
global IBM_IDCE IBM_GS IBM_INBOUNDCC IBM_MAX_WSTRIANG_SEG

ierr=1;


% BODINT_CELL:
for IG=1:N_GEOMETRY

   % The IG wet surface edges will be used to obtain intersections with grid planes on
   % increasing svar order.
   % Allocate BODINT_CELL:
   BODINT_CELL.NWSEGS = GEOM(IG).N_EDGES;
   BODINT_CELL.NWCROSS= zeros(1,BODINT_CELL.NWSEGS);
   BODINT_CELL.SVAR   = zeros(IBM_DELTA_NBCROSS,BODINT_CELL.NWSEGS);

   for IWSEDG=1:BODINT_CELL.NWSEGS

      % Seg Nodes location:
      SEG(NOD1:NOD2) = GEOM(IG).EDGES(NOD1:NOD2,IWSEDG);

      XYZ1(IAXIS:KAXIS) = GEOM(IG).VERTS(MAX_DIM*(SEG(NOD1)-1)+1:MAX_DIM*SEG(NOD1));
      XYZ2(IAXIS:KAXIS) = GEOM(IG).VERTS(MAX_DIM*(SEG(NOD2)-1)+1:MAX_DIM*SEG(NOD2));

      % Test if Segment lays on plane, If so drop, it was taken care of
      % previously: This is expensive think of switching to pointer X1FACEP.
      DROPSEG = false;
      for X1AXIS=IAXIS:KAXIS
         switch(X1AXIS)
           case(IAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 1.,  0., 0. ];
              X1FACE =  XFACE;
              DX1FACE = DXFACE;
              X1LO   = ILO_FACE-CCGUARD; X1HI   = IHI_FACE+CCGUARD;
           case(JAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 0.,  1., 0. ];
              X1FACE =  YFACE;
              DX1FACE = DYFACE;
              X1LO   = JLO_FACE-CCGUARD; X1HI   = JHI_FACE+CCGUARD;
           case(KAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 0.,  0., 1. ];
              X1FACE = ZFACE;
              DX1FACE = DZFACE;
              X1LO   = KLO_FACE-CCGUARD; X1HI   = KHI_FACE+CCGUARD;
         end

         % Optimized for UG:
         X1NOC=TRANS(NM).NOC(X1AXIS);
         MINX = min(XYZ1(X1AXIS),XYZ2(X1AXIS));
         MAXX = max(XYZ1(X1AXIS),XYZ2(X1AXIS));

         if (MAXX-MINX < GEOMEPS) % SEGMENT ALIGNED WITH PLANE.
            LSTR = X1LO; LEND = X1HI;
            if (X1NOC==0) % Optimized for Uniform Grid.
               LSTR = max(X1LO,  floor((MINX-GEOMEPS-X1FACE(X1LO))/DX1FACE(X1LO)) + X1LO);
               LEND = min(X1HI,   ceil((MAXX+GEOMEPS-X1FACE(X1LO))/DX1FACE(X1LO)) + X1LO);
            end
            for IPLN=LSTR:LEND
               X1PLN = X1FACE(IPLN);
               DROPSEG=(    abs(X1PLN-MAXX) < GEOMEPS) || ...
                       ((X1FACE(X1LO)-MAXX) > GEOMEPS) || ...
                       ((MAXX-X1FACE(X1HI)) > GEOMEPS);
               if (DROPSEG); break; end
            end
         end
         if (DROPSEG)
            break % break X1AXIS=IAXIS:KAXIS LOOP
         end
      end
      if (DROPSEG); continue; end

      % Edge length and tangent unit vector:
      DV(IAXIS:KAXIS) = XYZ2(IAXIS:KAXIS) - XYZ1(IAXIS:KAXIS);
      SLEN = sqrt( DV(IAXIS)^2. + DV(JAXIS)^2. + DV(KAXIS)^2. ); % Segment length.
      STANI(IAXIS:KAXIS) = DV(IAXIS:KAXIS) * SLEN^(-1.);         % Segment tangent versor.

      % Add segment ends as intersections:
      BODINT_CELL.NWCROSS(IWSEDG)  =    2; % Nodes 1,2 of segment are considered intersection.
      BODINT_CELL.SVAR(1,IWSEDG)   =    0; % Coordinate along stani of Node 1.
      BODINT_CELL.SVAR(2,IWSEDG)   = SLEN; % Coordinate along stani of Node 2.

      % Now find intersections:
      for X1AXIS=IAXIS:KAXIS
         switch(X1AXIS)
           case(IAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 1.,  0., 0. ];
              X1FACE =  XFACE;
              DX1FACE = DXFACE;
              X1LO   = ILO_FACE-CCGUARD; X1HI   = IHI_FACE+CCGUARD;
           case(JAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 0.,  1., 0. ];
              X1FACE =  YFACE;
              DX1FACE = DYFACE;
              X1LO   = JLO_FACE-CCGUARD; X1HI   = JHI_FACE+CCGUARD;
           case(KAXIS)
              PLNORMAL(IAXIS:KAXIS) = [ 0.,  0., 1. ];
              X1FACE = ZFACE;
              DX1FACE = DZFACE;
              X1LO   = KLO_FACE-CCGUARD; X1HI   = KHI_FACE+CCGUARD;
         end

         % Optimized for UG:
         X1NOC=TRANS(NM).NOC(X1AXIS);
         MINX = min(XYZ1(X1AXIS),XYZ2(X1AXIS));
         MAXX = max(XYZ1(X1AXIS),XYZ2(X1AXIS));
         LSTR = X1LO; LEND = X1HI;
         if (X1NOC==0) % Optimized for Uniform Grid.
            LSTR = max(X1LO,  floor((MINX-GEOMEPS-X1FACE(X1LO))/DX1FACE(X1LO)) + X1LO);
            LEND = min(X1HI,   ceil((MAXX+GEOMEPS-X1FACE(X1LO))/DX1FACE(X1LO)) + X1LO);
         end

         for IPLN=LSTR:LEND
            X1PLN = X1FACE(IPLN);
            OUTPLANE = ((X1PLN-MAXX) > GEOMEPS) || ((MINX-X1PLN) > GEOMEPS);
            if (OUTPLANE); continue; end % Make sure to drop jstr, jend if out of segment length.

            % Drop intersections in segment nodes:
            % Compute: dot(plnormal, xyzv - xypl):
            DOT1 = XYZ1(X1AXIS) - X1PLN;
            DOT2 = XYZ2(X1AXIS) - X1PLN;
            if (abs(DOT1) <= GEOMEPS); continue; end
            if (abs(DOT2) <= GEOMEPS); continue; end

            % Now regular case: Find svar and insert in BODINT_CELL.SVAR(:,IWSEDG):
            DENOM  = STANI(X1AXIS); % dot(stani,plnormal)
            PLANEEQ= DOT1   ;       % dot(xyz1(IAXIS:KAXIS),plnormal) - x1pln
            SVARI  = - PLANEEQ / DENOM;

            % Insertion sort, discard repeated intersections:
            NWCROSS = BODINT_CELL.NWCROSS(IWSEDG) + 1;
            BODINT_CELL.SVAR(NWCROSS,IWSEDG) = 1. / GEOMEPS;
            SAMEINT = false;
            for IBCR=1:NWCROSS
               if (abs(SVARI - BODINT_CELL.SVAR(IBCR,IWSEDG)) < GEOMEPS)
                  SAMEINT = true;
                  break
               end
               if ( SVARI  < BODINT_CELL.SVAR(IBCR,IWSEDG) ); break; end
            end
            if (SAMEINT); continue; end

            % Here copy from the back (updated nbcross) to the ibcr location:
            for IDUM = NWCROSS:-1:IBCR+1
               BODINT_CELL.SVAR(IDUM,IWSEDG) = BODINT_CELL.SVAR(IDUM-1,IWSEDG);
            end
            BODINT_CELL.SVAR(IBCR,IWSEDG) =   SVARI;
            BODINT_CELL.NWCROSS(IWSEDG)   = NWCROSS;
         end
      end
      
      % 3. The increasing svar intersections are used to define the INBOUNDCC type
      % cut-edges and Cartesian Cells containing them. Add to CUT_EDGE, define the
      % CUT_EDGE entry in CCVAR(i,j,k,IBM_IDCE):
      for IEDGE=1:BODINT_CELL.NWCROSS(IWSEDG)-1

         % Location along Segment:
         SVAR1 = BODINT_CELL.SVAR(IEDGE  ,IWSEDG);
         SVAR2 = BODINT_CELL.SVAR(IEDGE+1,IWSEDG);

         % Location of midpoint of cut-edge:
         SVAR12 = 0.5 * (SVAR1+SVAR2);

         % Define Cartesian cell this cut-edge belongs:
         % Optimized for UG version:
         XPOS = XYZ1(IAXIS) + SVAR12*STANI(IAXIS);
         if (TRANS(NM).NOC(IAXIS)==0)
            II2  = floor( (XPOS-XFACE(ILO_FACE))/DXFACE(ILO_FACE) ) + ILO_CELL;
            % Discard cut-edges on faces laying on x2hi and x3hi.
            if ( (II2 < ILO_CELL-CCGUARD) || (II2 > IHI_CELL+CCGUARD) ); continue; end
         else
            FOUND_SEG=false;
            for II2=ILO_CELL-CCGUARD:IHI_CELL+CCGUARD
               if ((XPOS-XFACE(II2-1)) >= 0. && (XFACE(II2)-XPOS) > 0.)
                  FOUND_SEG=true;
                  break
               end
            end
            if (~FOUND_SEG); continue; end
         end

         XPOS = XYZ1(JAXIS) + SVAR12*STANI(JAXIS);
         if (TRANS(NM).NOC(JAXIS)==0)
            JJ2  = floor( (XPOS-YFACE(JLO_FACE))/DYFACE(JLO_FACE) ) + JLO_CELL;
            if ( (JJ2 < JLO_CELL-CCGUARD) || (JJ2 > JHI_CELL+CCGUARD) ); continue; end
         else
            FOUND_SEG=false;
            for JJ2=JLO_CELL-CCGUARD:JHI_CELL+CCGUARD
               if ((XPOS-YFACE(JJ2-1)) >= 0. && (YFACE(JJ2)-XPOS) > 0.)
                  FOUND_SEG=true;
                  break
               end
            end
            if (~FOUND_SEG); continue; end
         end

         XPOS = XYZ1(KAXIS) + SVAR12*STANI(KAXIS);
         if (TRANS(NM).NOC(KAXIS)==0)
            KK2  = floor( (XPOS-ZFACE(KLO_FACE))/DZFACE(KLO_FACE) ) + KLO_CELL;
            if ( (KK2 < KLO_CELL-CCGUARD) || (KK2 > KHI_CELL+CCGUARD) ); continue; end
         else
            FOUND_SEG=false;
            for KK2=KLO_CELL-CCGUARD:KHI_CELL+CCGUARD
               if ((XPOS-ZFACE(KK2-1)) >= 0. && (ZFACE(KK2)-XPOS) > 0.)
                  FOUND_SEG=true;
                  break
               end
            end
            if (~FOUND_SEG); continue; end
         end


         % CCVAR edge number:
         if ( MESHES(NM).CCVAR(II2,JJ2,KK2,IBM_IDCE) > 0 ) % There is already
                                                                % an entry in CUT_EDGE.
            CEI = MESHES(NM).CCVAR(II2,JJ2,KK2,IBM_IDCE);
         else % We need a new entry in CUT_EDGE
            CEI = MESHES(NM).N_CUTEDGE_MESH + 1;
            MESHES(NM).N_CUTEDGE_MESH           = CEI;
            MESHES(NM).CCVAR(II2,JJ2,KK2,IBM_IDCE) = CEI;
            MESHES(NM).CUT_EDGE(CEI).NVERT  = 0;
            MESHES(NM).CUT_EDGE(CEI).NEDGE  = 0;
            MESHES(NM).CUT_EDGE(CEI).IJK(1:MAX_DIM+2) = [ II2, JJ2, KK2, 0, IBM_GS ];
            MESHES(NM).CUT_EDGE(CEI).STATUS = IBM_INBOUNDCC;
         end

         % Add vertices, non repeated vertex entries at this point.
         NVERT = MESHES(NM).CUT_EDGE(CEI).NVERT;
         % Define vertices for this segment:
         %               xv1                 yv1                 zv1
         XYZV1(IAXIS:KAXIS) = [ XYZ1(IAXIS)+SVAR1*STANI(IAXIS), ...
                                XYZ1(JAXIS)+SVAR1*STANI(JAXIS), ...
                                XYZ1(KAXIS)+SVAR1*STANI(KAXIS) ];
         [NVERT,INOD1]=INSERT_FACE_VERT(XYZV1,NM,CEI,NVERT);
         %               xv2                 yv2                 zv2
         XYZV2(IAXIS:KAXIS) = [ XYZ1(IAXIS)+SVAR2*STANI(IAXIS), ...
                                XYZ1(JAXIS)+SVAR2*STANI(JAXIS), ...
                                XYZ1(KAXIS)+SVAR2*STANI(KAXIS) ];
         [NVERT,INOD2]=INSERT_FACE_VERT(XYZV2,NM,CEI,NVERT);

         NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE + 1;
         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE) = [ INOD1, INOD2 ];

         MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,NEDGE) = [ INOD1, INOD2 ];
         MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,NEDGE) = ...
                                      [ GEOM(IG).EDGE_FACES(1,IWSEDG), ...
                                        GEOM(IG).EDGE_FACES(2,IWSEDG), ...
                                        GEOM(IG).EDGE_FACES(4,IWSEDG), IG ];
         MESHES(NM).CUT_EDGE(CEI).NVERT = NVERT;
         MESHES(NM).CUT_EDGE(CEI).NEDGE = NEDGE;
      end

   end
end

ierr=0;

return