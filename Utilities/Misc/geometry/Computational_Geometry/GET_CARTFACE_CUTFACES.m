function [ierr]=GET_CARTFACE_CUTFACES(NM,ISTR,IEND,JSTR,JEND,KSTR,KEND,BNDINT_FLAG)

global IAXIS JAXIS KAXIS CCGUARD IBM_SOLID IBM_CUTCFE NOD1 NOD2 NOD3 FCELL
global ILO_CELL IHI_CELL ILO_FACE IHI_FACE
global JLO_CELL JHI_CELL JLO_FACE JHI_FACE
global KLO_CELL KHI_CELL KLO_FACE KHI_FACE
global XFACE YFACE ZFACE
global MESHES IBM_UNDEFINED GEOMEPS MAX_DIM IBM_GASPHASE
global IBM_EGSC IBM_IDCE IBM_FGSC IBM_IDCF
global IBM_MAXCEELEM_FACE IBM_MAXVERTS_FACE
global IJK_COUNTED

global SEG_FACE_F SEG_FACE2_F XYZVERT_F NFACE_F CFELEM_F NVERT_F ANGSEG_F

ierr=1;


% Main Loop on block NM:
if (BNDINT_FLAG)
   IJK_COUNTED=zeros(IEND,JEND,KEND,KAXIS);
   BNDINT_LOW  = 1;
   BNDINT_HIGH = 3;
else
   if (CCGUARD==0); return; end
   BNDINT_LOW  = 4;
   BNDINT_HIGH = 4;
end

for IBNDINT=BNDINT_LOW:BNDINT_HIGH % 1,2 refers to block boundary faces, 3 to internal faces,
                                   % 4 guard-cell faces.

   % When switching to internal faces, copy number of external faces already computed.
   if (IBNDINT == 3)
       MESHES(NM).N_BBCUTFACE_MESH = MESHES(NM).N_CUTFACE_MESH;
   end

   for X1AXIS=IAXIS:KAXIS

      switch(X1AXIS)
      case(IAXIS)
         X2AXIS = JAXIS;
         X3AXIS = KAXIS;
         % IAXIS gasphase cut-faces:
         JLO = JLO_CELL; JHI = JHI_CELL;
         KLO = KLO_CELL; KHI = KHI_CELL;
         switch(IBNDINT)
         case(1)
            ILO = ILO_FACE; IHI = ILO_FACE;
         case(2)
            ILO = IHI_FACE; IHI = IHI_FACE;
         case(3)
            ILO = ILO_FACE+1; IHI = IHI_FACE-1;
         case(4)
            ILO = ILO_FACE-CCGUARD; IHI= IHI_FACE+CCGUARD;
            JLO = JLO-CCGUARD; JHI = JHI+CCGUARD;
            KLO = KLO-CCGUARD; KHI = KHI+CCGUARD;
         end
         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS;
         % Local indexing in x1, x2, x3:
         X1LO = ILO; X1HI = IHI;
         X2LO = JLO; X2HI = JHI;
         X3LO = KLO; X3HI = KHI;
         % Face coordinates in x1,x2,x3 axes:
         X1FACE = XFACE;
         X2FACE = YFACE;
         X3FACE = ZFACE;
      case(JAXIS)
         X2AXIS = KAXIS;
         X3AXIS = IAXIS;
         % JAXIS gasphase cut-faces:
         ILO = ILO_CELL; IHI = IHI_CELL;
         KLO = KLO_CELL; KHI = KHI_CELL;
         switch(IBNDINT)
         case(1)
            JLO = JLO_FACE; JHI = JLO_FACE;
         case(2)
            JLO = JHI_FACE; JHI = JHI_FACE;
         case(3)
            JLO = JLO_FACE+1; JHI = JHI_FACE-1;
         case(4)
            JLO = JLO_FACE-CCGUARD; JHI = JHI_FACE+CCGUARD;
            ILO = ILO-CCGUARD; IHI = IHI+CCGUARD;
            KLO = KLO-CCGUARD; KHI = KHI+CCGUARD;
         end
         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS;
         % Local indexing in x1, x2, x3:
         X1LO = JLO; X1HI = JHI;
         X2LO = KLO; X2HI = KHI;
         X3LO = ILO; X3HI = IHI;
         % Face coordinates in x1,x2,x3 axes:
         X1FACE = YFACE;
         X2FACE = ZFACE;
         X3FACE = XFACE;
      case(KAXIS)
         X2AXIS = IAXIS;
         X3AXIS = JAXIS;
         % KAXIS gasphase cut-faces:
         ILO = ILO_CELL; IHI = IHI_CELL;
         JLO = JLO_CELL; JHI = JHI_CELL;
         switch(IBNDINT)
         case(1)
            KLO = KLO_FACE; KHI = KLO_FACE;
         case(2)
            KLO = KHI_FACE; KHI = KHI_FACE;
         case(3)
            KLO = KLO_FACE+1; KHI = KHI_FACE-1;
         case(4)
            KLO = KLO_FACE-CCGUARD; KHI = KHI_FACE+CCGUARD;
            ILO = ILO-CCGUARD; IHI = IHI+CCGUARD;
            JLO = JLO-CCGUARD; JHI = JHI+CCGUARD;
         end
         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS;
         % Local indexing in x1, x2, x3:
         X1LO = KLO; X1HI = KHI;
         X2LO = ILO; X2HI = IHI;
         X3LO = JLO; X3HI = JHI;
         % Face coordinates in x1,x2,x3 axes:
         X1FACE = ZFACE;
         X2FACE = XFACE;
         X3FACE = YFACE;
      end

      % Loop on Cartesian faces, local x1, x2, x3 indexes:
      for II=X1LO:X1HI
         for KK=X3LO:X3HI
            for JJ=X2LO:X2HI

             % Face indexes:
             INDXI(IAXIS:KAXIS) = [ II, JJ, KK ]; % Local x1,x2,x3
             INDI = INDXI(XIAXIS);
             INDJ = INDXI(XJAXIS);
             INDK = INDXI(XKAXIS);

             % Drop if cut-face has already been counted:
             if( IJK_COUNTED(INDI,INDJ,INDK,X1AXIS) ); continue; end
             IJK_COUNTED(INDI,INDJ,INDK,X1AXIS)=true;
             if(MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_FGSC,X1AXIS) == IBM_SOLID); continue; end

             % Drop if face not cut-face:
             % Test for FACE Cartesian edges being cut:
             % If outface1 is true -> All regular edges for this face:
             % Edge at index KK-FCELL:
             INDXI1(IAXIS:KAXIS) = [ II, JJ, KK-FCELL ]; % Local x1,x2,x3
             INDI1 = INDXI1(XIAXIS);
             INDJ1 = INDXI1(XJAXIS);
             INDK1 = INDXI1(XKAXIS);
             % Edge at index KK-FCELL+1:
             INDXI2(IAXIS:KAXIS) = [ II, JJ, KK-FCELL+1 ]; % Local x1,x2,x3
             INDI2 = INDXI2(XIAXIS);
             INDJ2 = INDXI2(XJAXIS);
             INDK2 = INDXI2(XKAXIS);
             % Edge at index JJ-FCELL:
             INDXI3(IAXIS:KAXIS) = [ II, JJ-FCELL, KK ]; % Local x1,x2,x3
             INDI3 = INDXI3(XIAXIS);
             INDJ3 = INDXI3(XJAXIS);
             INDK3 = INDXI3(XKAXIS);
             % Edge at index jj-FCELL+1:
             INDXI4(IAXIS:KAXIS) = [ II, JJ-FCELL+1, KK ]; % Local x1,x2,x3
             INDI4 = INDXI4(XIAXIS);
             INDJ4 = INDXI4(XJAXIS);
             INDK4 = INDXI4(XKAXIS);

             OUTFACE1 = (MESHES(NM).ECVAR(INDI1,INDJ1,INDK1,IBM_EGSC,X2AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI2,INDJ2,INDK2,IBM_EGSC,X2AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI3,INDJ3,INDK3,IBM_EGSC,X3AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI4,INDJ4,INDK4,IBM_EGSC,X3AXIS) ~= IBM_CUTCFE);

             % Test for face with INB edges:
             % If outface2 is true -> no INB Edges associated with this face:
             OUTFACE2 = (MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_IDCE,X1AXIS) <= 0);

             % Drop if outface1 & outface2
             if (OUTFACE1 && OUTFACE2)
                % Test if IBM_FSID is SOLID:
                if ((MESHES(NM).ECVAR(INDI1,INDJ1,INDK1,IBM_EGSC,X2AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI2,INDJ2,INDK2,IBM_EGSC,X2AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI3,INDJ3,INDK3,IBM_EGSC,X3AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI4,INDJ4,INDK4,IBM_EGSC,X3AXIS) == IBM_SOLID) )
                   MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_FGSC,X1AXIS) = IBM_SOLID;
                end
                continue
             end

             MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_FGSC,X1AXIS)   = IBM_CUTCFE;

             % Build segment list:
             NSEG      = 0;
             NVERT     = 0;
             NFACE     = 0;

             SEG_FACE (NOD1:NOD2,1:IBM_MAXCEELEM_FACE)             = IBM_UNDEFINED;
             XYZVERT(IAXIS:KAXIS,1:IBM_MAXVERTS_FACE)              = 0.;
             ANGSEG(1:IBM_MAXCEELEM_FACE)                          = 0.;
             BODNUM(1:IBM_MAXCEELEM_FACE)                          = 1000000000;
             SEGTYPE(1:IBM_MAXCEELEM_FACE)                         = 0;

             % 1. Cartesian IBM_GASPHASE edges, cut-edges:
             % a. Make a list of segments:
             % Low x2 cut-edges:
             INDLC(IAXIS:KAXIS) = INDXI3(IAXIS:KAXIS);
             IEDG=INDI3; JEDG=INDJ3; KEDG=INDK3;
             CEI = MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_IDCE,X3AXIS);
             if ( CEI == 0 ) % Regular Edge, build segment from grid:
                if (MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_EGSC,X3AXIS) ~= IBM_SOLID)
                   % x,y,z of node 1:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)), ...
                                          X2FACE(INDLC(JAXIS)), ...
                                          X3FACE(INDLC(KAXIS)-FCELL+1) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)), ...
                                          X2FACE(INDLC(JAXIS)), ...
                                          X3FACE(INDLC(KAXIS)-FCELL) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = - pi / 2.;
                end
             else % Cut-edge, load CUT_EDGE(CEI) segments
                NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                for IEDGE=1:NEDGE
                   SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE);

                   % x,y,z of node 1:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2));
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1));
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = - pi / 2.;
                end
             end

             % High x2 cut-edges:
             INDLC(IAXIS:KAXIS) = INDXI4(IAXIS:KAXIS);
             IEDG=INDI4; JEDG=INDJ4; KEDG=INDK4;
             CEI = MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_IDCE,X3AXIS);
             if ( CEI == 0 ) % Regular Edge, build segment from grid:
                if (MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_EGSC,X3AXIS) ~= IBM_SOLID)
                   % x,y,z of node 1:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)), ...
                                          X2FACE(INDLC(JAXIS)), ...
                                          X3FACE(INDLC(KAXIS)-FCELL) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)), ...
                                          X2FACE(INDLC(JAXIS)), ...
                                          X3FACE(INDLC(KAXIS)-FCELL+1) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) =   pi / 2.;
                end
             else % Cut-edge, load CUT_EDGE(CEI) segments
                NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                for IEDGE=1:NEDGE
                   SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE);

                   % x,y,z of node 1:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1));
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2));
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) =   pi / 2.;
                end
             end

             % Low  x3 cut-edges:
             INDLC(IAXIS:KAXIS) = INDXI1(IAXIS:KAXIS);
             IEDG=INDI1; JEDG=INDJ1; KEDG=INDK1;
             CEI = MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_IDCE,X2AXIS);
             if ( CEI == 0 ) % Regular Edge, build segment from grid:
                if (MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_EGSC,X2AXIS) ~= IBM_SOLID)
                   % x,y,z of node 1:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)),       ...
                                          X2FACE(INDLC(JAXIS)-FCELL), ...
                                          X3FACE(INDLC(KAXIS)) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)),         ...
                                          X2FACE(INDLC(JAXIS)-FCELL+1), ...
                                          X3FACE(INDLC(KAXIS)) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = 0.;
                end
             else % Cut-edge, load CUT_EDGE(CEI) segments
                NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                for IEDGE=1:NEDGE
                   SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE);

                   % x,y,z of node 1:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1));
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2));
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = 0.;
                end
             end

             % High x3 cut-edges:
             INDLC(IAXIS:KAXIS) = INDXI2(IAXIS:KAXIS);
             IEDG=INDI2; JEDG=INDJ2; KEDG=INDK2;
             CEI = MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_IDCE,X2AXIS);
             if ( CEI == 0 ) % Regular Edge, build segment from grid:
                if (MESHES(NM).ECVAR(IEDG,JEDG,KEDG,IBM_EGSC,X2AXIS) ~= IBM_SOLID)
                   % x,y,z of node 1:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)),          ...
                                           X2FACE(INDLC(JAXIS)-FCELL+1), ...
                                           X3FACE(INDLC(KAXIS)) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZLC(IAXIS:KAXIS) = [ X1FACE(INDLC(IAXIS)),       ...
                                          X2FACE(INDLC(JAXIS)-FCELL), ...
                                          X3FACE(INDLC(KAXIS)) ];
                   X1 = XYZLC(XIAXIS);
                   X2 = XYZLC(XJAXIS);
                   X3 = XYZLC(XKAXIS);
                   XYZV(IAXIS:KAXIS) = [ X1, X2, X3 ];
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = pi;
                end
             else % Cut-edge, load CUT_EDGE(CEI) segments
                NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                for IEDGE=1:NEDGE
                   SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE);

                   % x,y,z of node 1:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2));
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1));
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   NSEG = NSEG + 1;
                   SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                   ANGSEG(NSEG) = pi;
                end
             end

             % Store Segment and Vertex list from Cartesian face boundary:
             XYZVERT_CART(IAXIS:KAXIS,1:NVERT)= XYZVERT(IAXIS:KAXIS,1:NVERT);
             SEG_FACE_CART(NOD1:NOD2,1:NSEG)  = SEG_FACE(NOD1:NOD2,1:NSEG);
             NVERT_CART=NVERT; NSEG_CART = NSEG;

             % 2. IBM_INBOUNDARY cut-edges assigned to this face:
             CEI = MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_IDCE,X1AXIS);
             if ( CEI > 0 ) % There are inboundary cut-edges
                NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;
                for IEDGE=1:NEDGE
                   SEG(NOD1:NOD2) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,IEDGE);
                   IBOD           = MESHES(NM).CUT_EDGE(CEI).INDSEG(4,IEDGE);
                   STYPE          = MESHES(NM).CUT_EDGE(CEI).INDSEG(5,IEDGE);

                   % x,y,z of node 1:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD2));
                   [NVERT,INOD1,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % x,y,z of node 2:
                   XYZV(IAXIS:KAXIS) = MESHES(NM).CUT_EDGE(CEI).XYZVERT(IAXIS:KAXIS,SEG(NOD1));
                   [NVERT,INOD2,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_FACE,XYZV,NVERT,XYZVERT);

                   % ADD segment:
                   VEC(NOD1:NOD2) = [ INOD1, INOD2 ];

                   % Insertion ADD segment:
                   INLIST = false;
                   for IDUM = 1:NSEG
                       if( (SEG_FACE(NOD1,IDUM)==VEC(NOD1) && SEG_FACE(NOD2,IDUM)==VEC(NOD2)) )
                           if (STYPE >=SEGTYPE(IDUM) && BODNUM(IDUM) > IBOD)
                               BODNUM(IDUM)=IBOD;
                               SEGTYPE(IDUM)=STYPE;
                           end
                           INLIST = true;
                           break
                       end
                       if( (SEG_FACE(NOD2,IDUM)==VEC(NOD1) && SEG_FACE(NOD1,IDUM)==VEC(NOD2)) )
                           if (STYPE >=SEGTYPE(IDUM) && BODNUM(IDUM) > IBOD)
                               SEG_FACE(NOD1:NOD2,IDUM) = VEC(NOD1:NOD2)';
                               BODNUM(IDUM)  =IBOD;
                               SEGTYPE(IDUM) =STYPE;
                           end
                           INLIST = true;
                           break
                       end
                   end
                   if (~INLIST)
                       NSEG = NSEG + 1;
                       SEG_FACE(NOD1:NOD2,NSEG) = [ INOD1, INOD2 ];
                       BODNUM(NSEG)             = IBOD;
                       SEGTYPE(NSEG)            = STYPE;
                       DX3 = XYZVERT(X3AXIS,INOD2)-XYZVERT(X3AXIS,INOD1);
                       DX2 = XYZVERT(X2AXIS,INOD2)-XYZVERT(X2AXIS,INOD1);
                       ANGSEG(NSEG) = atan2(DX3,DX2);

                   end

                end
             end

             % Here expand SEG_FACE to contain all halfedges of STYPE=1:
             COUNT=0;
             SEG_FACEAUX (NOD1:NOD2,1:IBM_MAXCEELEM_FACE)             = IBM_UNDEFINED;
             ANGSEGAUX(1:IBM_MAXCEELEM_FACE)                          = 0.;
             for ISEG=1:NSEG
                 COUNT=COUNT+1;
                 SEG_FACEAUX (NOD1:NOD2,COUNT) = SEG_FACE(NOD1:NOD2,ISEG);
                 ANGSEGAUX(COUNT) = ANGSEG(ISEG);
                 if(SEGTYPE(ISEG)==1)
                     COUNT=COUNT+1;
                     SEG_FACEAUX (NOD1:NOD2,COUNT) = SEG_FACE(NOD2:-1:NOD1,ISEG);
                     if(ANGSEG(ISEG) > 0.)
                        ANGSEGAUX(COUNT) = ANGSEG(ISEG)-pi;
                     else
                        ANGSEGAUX(COUNT) = ANGSEG(ISEG)+pi;
                     end
                 end
             end
             SEG_FACE = SEG_FACEAUX;
             ANGSEG   = ANGSEGAUX;
             NSEG=COUNT;

%              if(X1AXIS==KAXIS && INDI==11 && INDJ==18 && INDK==16)
%                  disp(['1 Cartface cutface:' num2str([INDI INDJ INDK X1AXIS])])
%                  figure
%                  subplot(1,3,1)
%                  hold on
%                  a=0.001;
%                  for ISEG=1:NSEG
%                      XYZSEG(:,NOD1)=XYZVERT(:,SEG_FACE(NOD1,ISEG));
%                      XYZSEG(:,NOD2)=XYZVERT(:,SEG_FACE(NOD2,ISEG));
%                      plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),'r','Marker','o')
%                  end
%                  for IVERT=1:NVERT
%                      text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                           XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
% 
%                  end
%                  xlabel('X')
%                  ylabel('Y')
%                  zlabel('Z')
%                  axis equal
%                  box on
%                  view([45 45])
%                  pause
% 
%              end
% 
             NOTDONE = true;
             while(NOTDONE)
                NOTDONE = false;
                % Counts edges that reach nodes:
                NUMEDG_NODE(1:IBM_MAXVERTS_FACE) = 0;
                for ISEG=1:NSEG
                   for II2=NOD1:NOD2
                      INOD = SEG_FACE(II2,ISEG);
                      NUMEDG_NODE(INOD) = NUMEDG_NODE(INOD) + 1;
                   end
                end

                % Drop segments with NUMEDG_NODE(INOD)=1:
                % The assumption here is that they are IBM_GG IBM_INBOUNDCF
                % segments with one node inside the Cartface i.e. case Fig
                % 9(a) in the CompGeom3D notes):
                COUNT = 0;
                SEG_FACEAUX (NOD1:NOD2,1:IBM_MAXCEELEM_FACE)             = IBM_UNDEFINED;
                ANGSEGAUX(1:IBM_MAXCEELEM_FACE)                          = 0.;
                for ISEG=1:NSEG
                   NUMNOD1 = NUMEDG_NODE(SEG_FACE(NOD1,ISEG));
                   NUMNOD2 = NUMEDG_NODE(SEG_FACE(NOD2,ISEG));
                   if ((NUMNOD1 > 1) && (NUMNOD2 > 1))
                      COUNT = COUNT + 1;
                      SEG_FACEAUX(NOD1:NOD2,COUNT) = SEG_FACE(NOD1:NOD2,ISEG);
                      ANGSEGAUX(COUNT) = ANGSEG(ISEG);
                   else
                      NOTDONE = true;
                   end
                end
                NSEG = COUNT;
                SEG_FACE = SEG_FACEAUX;
                ANGSEG   = ANGSEGAUX;
             end

             % Discard face with no conected edges or two edges connecting same nodes:
             if ( NSEG == 0 || (NSEG ==2 && ( any(SEG_FACE(NOD1:NOD2,1)==SEG_FACE(NOD2,2)) && ...
                                              any(SEG_FACE(NOD1:NOD2,1)==SEG_FACE(NOD1,2)) )) )
                MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_FGSC,X1AXIS) = IBM_SOLID;
                continue
             end

             % Add segments which have both ends attached to more than two segs:
             COUNT = 0;
             for ISEG=1:NSEG
                NUMNOD1 = NUMEDG_NODE(SEG_FACE(NOD1,ISEG));
                NUMNOD2 = NUMEDG_NODE(SEG_FACE(NOD2,ISEG));
                if ((NUMNOD1 > 2) && (NUMNOD2 > 2))
                   COUNT = COUNT + 1;
                   SEG_FACE(NOD1:NOD2,NSEG+COUNT) = SEG_FACE( [ NOD2, NOD1 ] ,ISEG);
                   if (ANGSEG(ISEG) >= 0.)
                      ANGSEG(NSEG+COUNT) = ANGSEG(ISEG) - pi;
                   else
                      ANGSEG(NSEG+COUNT) = ANGSEG(ISEG) + pi;
                   end
                end
             end
             NSEG = NSEG + COUNT;

             % Fill NODEDG_FACE(IEDGE,INOD), where iedge are edges
             % that contain inod as first node. This assumes edges are
             % ordered using the right hand rule on x2-x3 plane.
             % Also compute the edges angles in x2-x3 plane:
             % Reallocate NODEDG_FACE if NSEG+1 > SIZE_EDGES_NODEDG, or NVERT > SIZE_VERTS_NODEDG:
             %CALL REALLOCATE_NODEDG_FACE(NSEG,NVERT)
             NODEDG_FACE=zeros(NSEG+1,NVERT);
             for ISEG=1:NSEG
               INOD1 = SEG_FACE(NOD1,ISEG);
               NEDI  = NODEDG_FACE(1,INOD1) + 1; % Increase number of edges connected to node by 1.
               NODEDG_FACE(     1,INOD1) = NEDI;
               NODEDG_FACE(NEDI+1,INOD1) = ISEG;
             end

%              if(X1AXIS==KAXIS && INDI==11 && INDJ==18 && INDK==16)
%                  disp(['2 Cartface cutface:' num2str([INDI INDJ INDK X1AXIS])])
%                  subplot(1,3,2)
%                  hold on
%                  a=0.001;
%                  xlabel('X')
%                  ylabel('Y')
%                  zlabel('Z')
%                  axis equal
%                  box on
%                  view([45 45])
% 
%                  for ISEG=1:NSEG
%                      XYZSEG(:,NOD1)=XYZVERT(:,SEG_FACE(NOD1,ISEG));
%                      XYZSEG(:,NOD2)=XYZVERT(:,SEG_FACE(NOD2,ISEG));
%                      plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),'r','Marker','o')
%                      text(0.5*sum(XYZSEG(IAXIS,:))+a,0.5*sum(XYZSEG(JAXIS,:))+a,...
%                           0.5*sum(XYZSEG(KAXIS,:))+a,num2str(ANGSEG(ISEG)),'FontSize',10)
%                      pause
%                  end
%                  for IVERT=1:NVERT
%                      text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                           XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
% 
%                  end
%                  pause
% 
%              end


             % Now Reorder Segments, do tests:
             SEG_FACE2(NOD1:NOD3,1:IBM_MAXCEELEM_FACE) = IBM_UNDEFINED;  % [INOD1 INOD2 ICF]
             SEG_FLAG(1:IBM_MAXCEELEM_FACE) = true;

             ICF  = 1;
             ISEG = 1;
             NEWSEG = ISEG;
             COUNT= 1;
             CTSTART=COUNT;
             SEG_FACE2(NOD1:NOD3,COUNT) = [ SEG_FACE(NOD1,NEWSEG), SEG_FACE(NOD2,NEWSEG), ICF ];
             SEG_FLAG(ISEG) = false;
             NSEG_LEFT      = NSEG - 1;


             % Infamous infinite loop:
             while(1)

                FOUNDSEG = false;
                N2COUNT  = SEG_FACE2(NOD2,COUNT); % Node 2 of segment COUNT.
                ANGCOUNT = ANGSEG(NEWSEG);

                % Find Segment starting on Node 2 with smaller ANGSEG respect to COUNT.
                DANG = -1. / GEOMEPS;
                for ISS=2:NODEDG_FACE(1,N2COUNT)+1
                   ISEG = NODEDG_FACE(ISS,N2COUNT);
                   if ( SEG_FLAG(ISEG) ) % This seg hasn't been added to SEG_FACE2
                                         % Drop if seg is the opposite of count seg:
                      if ( SEG_FACE2(NOD1,COUNT) == SEG_FACE(NOD2,ISEG) && NUMEDG_NODE(N2COUNT) > 2); continue; end
                      DANGI = ANGSEG(ISEG) - ANGCOUNT;
                      if ( DANGI < 0. ); DANGI = DANGI + 2. * pi; end

                      if ( DANGI > DANG )
                         NEWSEG   =  ISEG;
                         DANG     = DANGI;
                         FOUNDSEG = true;
                      end
                   end
                end

                % Found a seg add to SEG_FACE2:
                if ( FOUNDSEG )
                   COUNT          = COUNT + 1;
                   SEG_FACE2(NOD1:NOD3,COUNT) = [ SEG_FACE(NOD1,NEWSEG), SEG_FACE(NOD2,NEWSEG), ICF ];
                   SEG_FLAG(NEWSEG) = false;
                   NSEG_LEFT      = NSEG_LEFT - 1;
                end

                % Test if line has closed on point shared any other cutface:
                if ( SEG_FACE2(NOD2,COUNT) == SEG_FACE2(NOD1,CTSTART) )
                   % Go for new cut-face on this Cartesian face.
                elseif ( FOUNDSEG )
                   continue
                end

                % Break loop:
                if ( NSEG_LEFT == 0 ); break; end

                % Start a new cut-face on this Cartesian face:
                ICF = ICF + 1;
                for ISEG=1:NSEG
                   if ( SEG_FLAG(ISEG) )
                      COUNT  = COUNT + 1;
                      CTSTART= COUNT;
                      SEG_FACE2(NOD1:NOD3,COUNT) = [ SEG_FACE(NOD1,ISEG), SEG_FACE(NOD2,ISEG), ICF ];
                      SEG_FLAG(ISEG) = false;
                      NSEG_LEFT      = NSEG_LEFT - 1;
                      break
                   end
                end

             end

             % Load ordered nodes to CFELEM:
             NFACE = ICF;
             % Reallocate CFELEM ARRAY if necessary:
             % CALL REALLOCATE_LOCAL_CFELEM(NSEG,NFACE)
             CFELEM(NSEG,NFACE) = IBM_UNDEFINED;
             for ICF=1:NFACE
                NP = 0;
                for ISEG=1:NSEG
                   if ( SEG_FACE2(NOD3,ISEG) == ICF )
                      NP = NP + 1;
                      CFELEM(1,ICF)    = NP;
                      CFELEM(NP+1,ICF) = SEG_FACE2(NOD1,ISEG);
                   end
                end
             end

             [rows,cols]=size(CFELEM);
             CFELEM2 = IBM_UNDEFINED*ones(rows,cols);
             NP=0;
             for ICF=1:NFACE
                if(CFELEM(1,ICF)>2)
                   NP=NP+1;
                   CFELEM2(:,NP) = CFELEM(:,ICF);
                end
             end
             CFELEM = CFELEM2;
             NFACE = NP;


%              if(X1AXIS==KAXIS && INDI==11 && INDJ==18 && INDK==16)
%                  disp(['3 Cartface cutface:' num2str([INDI INDJ INDK X1AXIS])])
%                  subplot(1,3,3)
%                  axis equal; box on;
%                  view([45 45])
%                  xlabel('X')
%                  ylabel('Y')
%                  zlabel('Z')
%                  hold on
%                  for JCF=1:NFACE
%                      NELEM   = CFELEM(1,JCF)
%                      CFELEM2 = CFELEM(2:NELEM+1,JCF);
%                      [hp]=patch(XYZVERT(IAXIS,CFELEM2),XYZVERT(JAXIS,CFELEM2),XYZVERT(KAXIS,CFELEM2),'b','Marker','o');
%                      set(hp,'FaceAlpha',0.3)
%                      pause
%                  end
%                  a=0.0005;
%                  for IVERT=1:NVERT
%                      text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                           XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
%                  end
%                  pause
% 
%              end
% 
%              if(X1AXIS==KAXIS && INDI==11 && INDJ==18 && INDK==16)
%                  disp(['Face=' num2str(MESHES(NM).N_CUTFACE_MESH + MESHES(NM).N_GCCUTFACE_MESH + 1)])
%                  NVERT_F = NVERT;
%                  NFACE_F = NFACE;
%                  SEG_FACE_F = SEG_FACE;
%                  ANGSEG_F   = ANGSEG;
%                  SEG_FACE2_F = SEG_FACE2;
%                  
%                  XYZVERT_F   = XYZVERT;
%                  CFELEM_F    = CFELEM;
%              end

             % Compute area and Centroid, in local x1, x2, x3 coords:
             AREAV(1:NFACE)                 = 0.;
             XYZCEN(IAXIS:KAXIS,1:NFACE)    = 0.;
             INXAREA(IAXIS:KAXIS,1:NFACE)   = 0.;
             INXSQAREA(IAXIS:KAXIS,1:NFACE) = 0.;
             for ICF=1:NFACE
                NP    = CFELEM(1,ICF);
                for IPT=2:NP+1
                   ICF_PT = CFELEM(IPT,ICF);
                   % Define closed Polygon centered in First Point:
                   XY(IAXIS:JAXIS,IPT-1) = [ XYZVERT(X2AXIS,ICF_PT)-XYZVERT(X2AXIS,CFELEM(2,ICF)), ...
                                             XYZVERT(X3AXIS,ICF_PT)-XYZVERT(X3AXIS,CFELEM(2,ICF)) ];
                end
                ICF_PT = CFELEM(2,ICF);
                XY(IAXIS:JAXIS,NP+1) = [ XYZVERT(X2AXIS,ICF_PT)-XYZVERT(X2AXIS,CFELEM(2,ICF)), ...
                                          XYZVERT(X3AXIS,ICF_PT)-XYZVERT(X3AXIS,CFELEM(2,ICF)) ];

                % Get Area and Centroid properties of Cut-face:
                AREA = 0.;
                for II2=1:NP
                   AREA = AREA + ( XY(IAXIS,II2) * XY(JAXIS,II2+1) - ...
                                   XY(JAXIS,II2) * XY(IAXIS,II2+1) );
                end
                AREA = AREA / 2.;
                % Now Centroids:
                % In x2:
                CX2 = 0.;
                for II2=1:NP
                   CX2 = CX2 + ( XY(IAXIS,II2)+XY(IAXIS,II2+1)) * ...
                               ( XY(IAXIS,II2)*XY(JAXIS,II2+1)  - ...
                                 XY(JAXIS,II2)*XY(IAXIS,II2+1) );
                end
                CX2 = CX2 / (6. * AREA) + XYZVERT(X2AXIS,CFELEM(2,ICF));
                % In x3:
                CX3 = 0.;
                for II2=1:NP
                   CX3 = CX3 + ( XY(JAXIS,II2)+XY(JAXIS,II2+1)) * ...
                               ( XY(IAXIS,II2)*XY(JAXIS,II2+1)  - ...
                                 XY(JAXIS,II2)*XY(IAXIS,II2+1) );
                end
                CX3 = CX3 / (6. * AREA) + XYZVERT(X3AXIS,CFELEM(2,ICF));

                % Add to cut-face:
                AREAV(ICF) = AREA;
                XYZCEN(IAXIS:KAXIS,ICF) = [  X1FACE(II), CX2, CX3 ];

                % Fields for cut-cell volume/centroid computation:
                % dot(e1,nc)*int(x1)dA, where x=x1face(ii) constant and nc=e1:
                INXAREA(IAXIS,ICF) = 1. * X1FACE(II) * AREA;
                INXAREA(JAXIS,ICF) = 0.;
                INXAREA(KAXIS,ICF) = 0.;
                % dot(e1,nc)*int(x1^2)dA, where x=x1face(ii) constant and nc=e1:
                INXSQAREA(IAXIS,ICF) = 1. * X1FACE(II)^2. * AREA;
                % dot(e2,nc)*int(x2^2)dA, where nc=e1 => dot(e2,nc)=0:
                INXSQAREA(JAXIS,ICF) = 0.;
                % dot(e3,nc)*int(x3^2)dA, where nc=e1 => dot(e3,nc)=0:
                INXSQAREA(KAXIS,ICF) = 0.;

             end

             % Figure out if a cut-face is completely inside any of the
             % others (that is, it is a hole on the GASPHASE):
             FINFACE =  zeros(1,NFACE);
             NFACE2  =  NFACE;
             for ICF1=1:NFACE2
                % Test that ICF1 has a negative area (case of holes)
                AREA1 = AREAV(ICF1);
                if ( AREA1 < -GEOMEPS )
                   for ICF2=1:NFACE2
                      % Drop if same face:
                      if ( ICF1 == ICF2 ); continue; end

                      % Centroid node for ICF1:
                      XYC1(1:2) = XYZCEN( [ JAXIS, KAXIS ] , ICF1 ); % [x2axis x3axis]

                      % Polygon nodes for ICF2:
                      NP2 = CFELEM(1,ICF2);
                      for IPT=2:NP2+1
                         ICF_PT = CFELEM(IPT,ICF2);
                         % Define closed Polygon:
                         XY(IAXIS:JAXIS,IPT-1) = [ XYZVERT(X2AXIS,ICF_PT), XYZVERT(X3AXIS,ICF_PT) ];
                      end

                      [PTSFLAG]=TEST_PT_INPOLY(NP2,XY,XYC1);

                      if ( PTSFLAG ) % Centroid of face 1 inside Face 2.

                         FINFACE(ICF1) = ICF2;
                         NFACE = NFACE - 1;

                         % Redefine areas in case of faces with holes:
                         AREA2 = AREAV(ICF2);

                         % Area with hole, AREA1 has negative sign:
                         AREAH = AREA2 + AREA1;

                         if(abs(AREAH) < GEOMEPS)
                             FINFACE(ICF2) = ICF1;
                             continue
                         end

                         % Centroid with hole:
                         XYC2(1:2) = XYZCEN( [ JAXIS, KAXIS ] , ICF2 );  % [x2axis x3axis]
                         XYH(1:2)  = (AREA1 * XYC1(1:2) + AREA2 * XYC2(1:2)) / AREAH;

                         % So ICF2 has the area with hole properties:
                         AREAV(ICF2) = AREAH;
                         XYZCEN(JAXIS,ICF2) = XYH(IAXIS);
                         XYZCEN(KAXIS,ICF2) = XYH(JAXIS);

                         % Other geom variables:
                         INXAREA(IAXIS:KAXIS,ICF2)  =  INXAREA(IAXIS:KAXIS,ICF2)+  INXAREA(IAXIS:KAXIS,ICF1);
                         INXSQAREA(IAXIS:KAXIS,ICF2)=INXSQAREA(IAXIS:KAXIS,ICF2)+INXSQAREA(IAXIS:KAXIS,ICF1);

                         break
                      end
                   end
                end
             end

             % Now enhance CFELEM for faces with holes nodes:
             for ICF1=1:NFACE2
                ICF2 = FINFACE(ICF1);
                if ( ICF2 > 0 ) % Allows for up to one hole per IBM_GASPHASE cut-face.
                   % Load points
                   NP1    = CFELEM(1,ICF1);
                   NP2    = CFELEM(1,ICF2);
                   NP     = (NP1+1) + (NP2+1);

                   % Here reallocate CFELEM, CFE, CFEL if NP > SIZE_VERTS_CFELEM:
                   % CALL REALLOCATE_LOCAL_VERT_CFELEM(NP+1)
                   CFE(1) = NP;

                   for II2=2:NP1+1
                      CFE(II2) = CFELEM(II2,ICF1);
                   end
                   II2 = (NP1+1) + 1;
                   CFE(II2) = CFELEM(2,ICF1);

                   % Load last point location:
                   ILOC = 2;
                   DIST12 = 1. / GEOMEPS;
                   XYC1(1:2) = [ XYZVERT(X2AXIS,CFE(II2)), XYZVERT(X3AXIS,CFE(II2)) ];
                   for COUNT=2:NP2+1
                      XYC2(1:2) = [ XYZVERT(X2AXIS,CFELEM(COUNT,ICF2)), XYZVERT(X3AXIS,CFELEM(COUNT,ICF2)) ];
                      D12 = sqrt( (XYC1(1)-XYC2(1))^2. + (XYC1(2)-XYC2(2))^2. );
                      if( D12 < DIST12 )
                         DIST12 = D12;
                         ILOC = COUNT;
                      end
                   end
                   if (ILOC > 2)
                      % Rebuild CFELEM(:,ICF2) such that the first point is ILOC:
                      CFEL(2:2+(NP2+1)-ILOC)    = CFELEM(ILOC:NP2+1,ICF2);
                      CFEL(3+(NP2+1)-ILOC:NP2+1)= CFELEM(2:ILOC-1 ,ICF2);
                      CFELEM(2:NP2+1 ,ICF2)     = CFEL(2:NP2+1);
                   end

                   COUNT = 1;
                   for II2=(NP1+1)+2:(NP1+1)+1+NP2
                      COUNT    = COUNT + 1;
                      CFE(II2) = CFELEM(COUNT,ICF2);
                   end
                   II2 = NP + 1;
                   CFE(II2) = CFELEM(2,ICF2);

                   % Copy CFE into CFELEM(1:np+1,icf2):
                   CFELEM(1:NP+1,ICF2) = CFE(1:NP+1);

                end
             end

             NVERTFACE = max(CFELEM(1,1:NFACE)) + 1;

             % This is a cut-face, allocate space:
             NCUTFACE = MESHES(NM).N_CUTFACE_MESH + MESHES(NM).N_GCCUTFACE_MESH + 1;
             if (BNDINT_FLAG)
                MESHES(NM).N_CUTFACE_MESH = NCUTFACE;
             else
                MESHES(NM).N_GCCUTFACE_MESH = MESHES(NM).N_GCCUTFACE_MESH + 1;
             end
             MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_IDCF,X1AXIS) = NCUTFACE;
             MESHES(NM).CUT_FACE(NCUTFACE).NVERT  = NVERT;
             MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = NFACE;
             MESHES(NM).CUT_FACE(NCUTFACE).IJK(1:MAX_DIM+1) = [ INDI, INDJ, INDK, X1AXIS ];
             MESHES(NM).CUT_FACE(NCUTFACE).STATUS = IBM_GASPHASE;
             MESHES(NM).CUT_FACE(NCUTFACE).XYZVERT(IAXIS:KAXIS,1:NVERT) = XYZVERT(IAXIS:KAXIS,1:NVERT);

             % Load Ordered nodes to CFELEM and geom properties:
             COUNT = 0;
             for ICF=1:NFACE2
                if ( FINFACE(ICF) > 0 ); continue; end % icf is a hole on another cut-face.
                COUNT = COUNT + 1;
                % Connectivity:
                MESHES(NM).CUT_FACE(NCUTFACE).CFELEM(1:NVERTFACE,COUNT) = ...
                                              CFELEM(1:NVERTFACE, ICF);
                % Geom Properties:
                MESHES(NM).CUT_FACE(NCUTFACE).AREA(COUNT) = AREAV(ICF);
                MESHES(NM).CUT_FACE(NCUTFACE).XYZCEN(IAXIS:KAXIS,COUNT) = ...
                                              XYZCEN( [ XIAXIS, XJAXIS, XKAXIS ] ,ICF);

                % Fields for cut-cell volume/centroid computation:
                % dot(i,nc)*int(x)dA, where nc=j => dot(i,nc)=0:
                MESHES(NM).CUT_FACE(NCUTFACE).INXAREA(COUNT)   =   INXAREA(XIAXIS,ICF);
                % dot(i,nc)*int(x^2)dA, where nc=j => dot(i,nc)=0:
                MESHES(NM).CUT_FACE(NCUTFACE).INXSQAREA(COUNT) = INXSQAREA(XIAXIS,ICF);
                % dot(j,nc)*int(y^2)dA, where y=yface(J) constant nc=j:
                MESHES(NM).CUT_FACE(NCUTFACE).JNYSQAREA(COUNT) = INXSQAREA(XJAXIS,ICF);
                % dot(k,nc)*int(z^2)dA, where nc=j => dot(k,nc)=0:
                MESHES(NM).CUT_FACE(NCUTFACE).KNZSQAREA(COUNT) = INXSQAREA(XKAXIS,ICF);
             end
             % Final number of cut-faces in the gas region of the face:
             NFACE = COUNT;
             MESHES(NM).CUT_FACE(NCUTFACE).NFACE  = NFACE;

            end % JJ
         end % KK
      end % II
   end

end

if (BNDINT_FLAG)
   % Here we mark faces on the guard-cell region for the computaiton of grid aligned INBOUNDARY faces
   % on CARTCELL_CUTFACES to work correctly:
   for X1AXIS=IAXIS:KAXIS

      switch(X1AXIS)
      case(IAXIS)

         X2AXIS = JAXIS;
         X3AXIS = KAXIS;

         % IAXIS gasphase cut-faces:
         ILO = ILO_FACE-CCGUARD; IHI = IHI_FACE+CCGUARD;
         JLO = JLO_CELL-CCGUARD; JHI = JHI_CELL+CCGUARD;
         KLO = KLO_CELL-CCGUARD; KHI = KHI_CELL+CCGUARD;

         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS;

         % Local indexing in x1, x2, x3:
         X1LO = ILO; X1HI = IHI;
         X2LO = JLO; X2HI = JHI;
         X3LO = KLO; X3HI = KHI;

      case(JAXIS)

         X2AXIS = KAXIS;
         X3AXIS = IAXIS;

         % JAXIS gasphase cut-faces:
         JLO = JLO_FACE-CCGUARD; JHI = JHI_FACE+CCGUARD;
         ILO = ILO_CELL-CCGUARD; IHI = IHI_CELL+CCGUARD;
         KLO = KLO_CELL-CCGUARD; KHI = KHI_CELL+CCGUARD;

         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS;

         % Local indexing in x1, x2, x3:
         X1LO = JLO; X1HI = JHI;
         X2LO = KLO; X2HI = KHI;
         X3LO = ILO; X3HI = IHI;

      case(KAXIS)

         X2AXIS = IAXIS;
         X3AXIS = JAXIS;

         % KAXIS gasphase cut-faces:
         KLO = KLO_FACE-CCGUARD; KHI = KHI_FACE+CCGUARD;
         ILO = ILO_CELL-CCGUARD; IHI = IHI_CELL+CCGUARD;
         JLO = JLO_CELL-CCGUARD; JHI = JHI_CELL+CCGUARD;

         % location in I,J,K od x2,x2,x3 axes:
         XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS;

         % Local indexing in x1, x2, x3:
         X1LO = KLO; X1HI = KHI;
         X2LO = ILO; X2HI = IHI;
         X3LO = JLO; X3HI = JHI;

      end

      % Loop on Cartesian faces, local x1, x2, x3 indexes:
      for II=X1LO:X1HI
         for KK=X3LO:X3HI
            for JJ=X2LO:X2HI

             % Face indexes:
             INDXI(IAXIS:KAXIS) = [ II, JJ, KK ]; % Local x1,x2,x3
             INDI = INDXI(XIAXIS);
             INDJ = INDXI(XJAXIS);
             INDK = INDXI(XKAXIS);

             % Drop if cut-face has already been counted:
             if( IJK_COUNTED(INDI,INDJ,INDK,X1AXIS) ); continue; end

             % Drop if face not cut-face:
             % Test for FACE Cartesian edges being cut:
             % If outface1 is true -> All regular edges for this face:
             % Edge at index KK-FCELL:
             INDXI1(IAXIS:KAXIS) = [ II, JJ, KK-FCELL ]; % Local x1,x2,x3
             INDI1 = INDXI1(XIAXIS);
             INDJ1 = INDXI1(XJAXIS);
             INDK1 = INDXI1(XKAXIS);
             % Edge at index KK-FCELL+1:
             INDXI2(IAXIS:KAXIS) = [ II, JJ, KK-FCELL+1 ]; % Local x1,x2,x3
             INDI2 = INDXI2(XIAXIS);
             INDJ2 = INDXI2(XJAXIS);
             INDK2 = INDXI2(XKAXIS);
             % Edge at index JJ-FCELL:
             INDXI3(IAXIS:KAXIS) = [ II, JJ-FCELL, KK ]; % Local x1,x2,x3
             INDI3 = INDXI3(XIAXIS);
             INDJ3 = INDXI3(XJAXIS);
             INDK3 = INDXI3(XKAXIS);
             % Edge at index jj-FCELL+1:
             INDXI4(IAXIS:KAXIS) = [ II, JJ-FCELL+1, KK ]; % Local x1,x2,x3
             INDI4 = INDXI4(XIAXIS);
             INDJ4 = INDXI4(XJAXIS);
             INDK4 = INDXI4(XKAXIS);

             OUTFACE1 = (MESHES(NM).ECVAR(INDI1,INDJ1,INDK1,IBM_EGSC,X2AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI2,INDJ2,INDK2,IBM_EGSC,X2AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI3,INDJ3,INDK3,IBM_EGSC,X3AXIS) ~= IBM_CUTCFE) && ...
                        (MESHES(NM).ECVAR(INDI4,INDJ4,INDK4,IBM_EGSC,X3AXIS) ~= IBM_CUTCFE);

             % Test for face with INB edges:
             % If outface2 is true -> no INB Edges associated with this face:
             OUTFACE2 = (MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_IDCE,X1AXIS) <= 0);

             % Drop if outface1 & outface2
             if (OUTFACE1 && OUTFACE2)
                % Test if IBM_FSID is SOLID:
                if ((MESHES(NM).ECVAR(INDI1,INDJ1,INDK1,IBM_EGSC,X2AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI2,INDJ2,INDK2,IBM_EGSC,X2AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI3,INDJ3,INDK3,IBM_EGSC,X3AXIS) == IBM_SOLID) && ...
                    (MESHES(NM).ECVAR(INDI4,INDJ4,INDK4,IBM_EGSC,X3AXIS) == IBM_SOLID) )
                   MESHES(NM).FCVAR(INDI,INDJ,INDK,IBM_FGSC,X1AXIS) = IBM_SOLID;
                end
                continue
             end

            end % JJ
         end % KK
      end % II

   end

else
   clear IJK_COUNTED

end

ierr=0;

return
