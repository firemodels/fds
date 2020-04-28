function [INB_FLG,NVERT,XYVERT,NEDGE,CEELEM,INDSEG]=GET_TRIANG_FACE_INT(X2AXIS,X3AXIS,FVERT,CEI,NM)

global IAXIS JAXIS DELTA_VERT IBM_MAXCEELEM_FACE IBM_MAX_WSTRIANG_SEG
global MESHES IBM_UNDEFINED GEOMEPS IBM_MAXVERTS_FACE GEOM
global NOD1 NOD2 NOD3 NOD4 EDG1 EDG2 EDG3
global BODINT_PLANE LOW_IND HIGH_IND


% Default return values:
INB_FLG = false;
NVERT   = 0;
NEDGE   = 0;
SIZE_X2X3VERT = DELTA_VERT;
X2X3VERT      =zeros(JAXIS,DELTA_VERT);
XYVERT        =zeros(JAXIS,DELTA_VERT);
CEELEM        = IBM_UNDEFINED*ones(NOD2,IBM_MAXCEELEM_FACE);
INDSEG        = IBM_UNDEFINED*ones(IBM_MAX_WSTRIANG_SEG+2,IBM_MAXCEELEM_FACE);
ISTRIEDGE     = ones(1,IBM_MAXCEELEM_FACE);

if ( CEI ~= 0 )
   NVERT = MESHES(NM).CUT_EDGE(CEI).NVERT;
   NEDGE = MESHES(NM).CUT_EDGE(CEI).NEDGE;

   if (NVERT > SIZE_X2X3VERT)
      SIZE_X2X3VERT = NVERT + DELTA_VERT;
      X2X3VERT=zeros(JAXIS,SIZE_X2X3VERT);
   end

   X2X3VERT(IAXIS,1:NVERT)   = MESHES(NM).CUT_EDGE(CEI).XYZVERT(X2AXIS,1:NVERT);
   X2X3VERT(JAXIS,1:NVERT)   = MESHES(NM).CUT_EDGE(CEI).XYZVERT(X3AXIS,1:NVERT);

   CEELEM(NOD1:NOD2,1:NEDGE) = MESHES(NM).CUT_EDGE(CEI).CEELEM(NOD1:NOD2,1:NEDGE);
   INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,1:NEDGE) = ...
   MESHES(NM).CUT_EDGE(CEI).INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,1:NEDGE);
   MESHES(NM).CUT_EDGE(CEI).NEDGE1=NEDGE;
end

% Quick discard test:
X2FMIN = min(FVERT(IAXIS,NOD1:NOD4)); X2FMAX = max(FVERT(IAXIS,NOD1:NOD4));
X3FMIN = min(FVERT(JAXIS,NOD1:NOD4)); X3FMAX = max(FVERT(JAXIS,NOD1:NOD4));

% Loop in-plane Surface Elements:
INTEST = false;
for ITRI=1:BODINT_PLANE.NTRIS
    % Elements nodes location, in x2-x3 coordinates:
    TRI(NOD1:NOD3) = BODINT_PLANE.TRIS(NOD1:NOD3,ITRI)';
    for INOD=NOD1:NOD3
       XYEL(IAXIS:JAXIS,INOD) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] ,TRI(INOD));
    end
    OUTX2= ((X2FMIN-max(XYEL(IAXIS,NOD1:NOD3))) > GEOMEPS) || ...
           ((min(XYEL(IAXIS,NOD1:NOD3))-X2FMAX) > GEOMEPS); % Triang out of Face in x2 dir
    OUTX3= ((X3FMIN-max(XYEL(JAXIS,NOD1:NOD3))) > GEOMEPS) || ...
           ((min(XYEL(JAXIS,NOD1:NOD3))-X3FMAX) > GEOMEPS); % Triang out of Face in x3 dir
    OUTFACE = OUTX2 || OUTX3;
    if (~OUTFACE)
        INTEST = true;
        break
    end
end
% Run on Triangle edges found:
for ISEG=1:BODINT_PLANE.NSEGS
   SEG(NOD1:NOD2) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG)';
   for INOD=NOD1:NOD2
      XYEL(IAXIS:JAXIS,INOD) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] ,SEG(INOD));
   end
   OUTX2= ((X2FMIN-max(XYEL(IAXIS,NOD1:NOD2))) > GEOMEPS) || ...
          ((min(XYEL(IAXIS,NOD1:NOD2))-X2FMAX) > GEOMEPS); % Segment out of Face in x2 dir
   OUTX3= ((X3FMIN-max(XYEL(JAXIS,NOD1:NOD2))) > GEOMEPS) || ...
          ((min(XYEL(JAXIS,NOD1:NOD2))-X3FMAX) > GEOMEPS); % Segment out of Face in x3 dir
   OUTFACE = OUTX2 || OUTX3;
   if (~OUTFACE)
       INTEST = true;
       break
   end
end

if (~INTEST); return; end

% Now if intest is true figure out if there are triangles-face intersection
% Polygons:
NFVERT = 4;
NTVERT = 3;
NSVERT = 2;

% First Vertices:
FVERT_IN_TRIANG=zeros(NFVERT,BODINT_PLANE.NTRIS);
TRIVERT_IN_FACE=zeros(NTVERT,BODINT_PLANE.NTRIS);

NINTP = NVERT;
% Loop in-plane Surface Elements:
for ITRI=1:BODINT_PLANE.NTRIS

   NINTP_TRI =             0;
   TRINODS   = IBM_UNDEFINED*ones(1,IBM_MAXVERTS_FACE);

   % Elements nodes location, in x2-x3 coordinates:
   TRI(NOD1:NOD3) = BODINT_PLANE.TRIS(NOD1:NOD3,ITRI)';
   for INOD=NOD1:NOD3
      XYEL(IAXIS:JAXIS,INOD) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] ,TRI(INOD));
   end

   % Cycle if Triangles BBOX not intersecting face:
   OUTX2= ((X2FMIN-max(XYEL(IAXIS,NOD1:NOD3))) > GEOMEPS) || ...
          ((min(XYEL(IAXIS,NOD1:NOD3))-X2FMAX) > GEOMEPS); % Triang out of Face in x2 dir
   OUTX3= ((X3FMIN-max(XYEL(JAXIS,NOD1:NOD3))) > GEOMEPS) || ...
          ((min(XYEL(JAXIS,NOD1:NOD3))-X3FMAX) > GEOMEPS); % Triang out of Face in x3 dir
   OUTFACE = OUTX2 || OUTX3;
   if (OUTFACE); continue; end

   if (BODINT_PLANE.X1NVEC(ITRI) < 0) % ROTATE NODE 2 AND 3 LOCATIONS
      DUMMY(IAXIS:JAXIS)     = XYEL(IAXIS:JAXIS,NOD2)';
      XYEL(IAXIS:JAXIS,NOD2) = XYEL(IAXIS:JAXIS,NOD3);
      XYEL(IAXIS:JAXIS,NOD3) = DUMMY(IAXIS:JAXIS)';

      TSEGS(NOD1:NOD2,EDG1) = BODINT_PLANE.TRIS( [ 2, 1 ] ,ITRI);
      TSEGS(NOD1:NOD2,EDG2) = BODINT_PLANE.TRIS( [ 3, 2 ] ,ITRI);
      TSEGS(NOD1:NOD2,EDG3) = BODINT_PLANE.TRIS( [ 1, 3 ] ,ITRI);
   else
      TSEGS(NOD1:NOD2,EDG1) = BODINT_PLANE.TRIS( [ 1, 2 ] ,ITRI);
      TSEGS(NOD1:NOD2,EDG2) = BODINT_PLANE.TRIS( [ 2, 3 ] ,ITRI);
      TSEGS(NOD1:NOD2,EDG3) = BODINT_PLANE.TRIS( [ 3, 1 ] ,ITRI);
   end

   
   % a. Test if Triangles vertices Lay on Faces area, including face boundary:
   for IPT=1:NTVERT
      OUTX2= ((X2FMIN-XYEL(IAXIS,IPT)) > GEOMEPS) || ...
             ((XYEL(IAXIS,IPT)-X2FMAX) > GEOMEPS);  % Triang out of Face in x2 dir
      OUTX3= ((X3FMIN-XYEL(JAXIS,IPT)) > GEOMEPS) || ...
             ((XYEL(JAXIS,IPT)-X3FMAX) > GEOMEPS);  % Triang out of Face in x3 dir
      OUTFACE = OUTX2 || OUTX3;

      if ( OUTFACE ); continue; end

      % Insertion add point to intersection list:
      XP(IAXIS:JAXIS) = XYEL(IAXIS:JAXIS,IPT)';
      [NINTP,SIZE_X2X3VERT,X2X3VERT,INOD]=INSERT_POINT_2D(XP,NINTP,SIZE_X2X3VERT,X2X3VERT);

      % Insert sort node to triangles local list
      TRUETHAT = true;
      for INP=1:NINTP_TRI
         if (TRINODS(INP) == INOD)
            TRUETHAT = false;
            break
         end
      end
      if ( TRUETHAT ) % new inod entry on list
         NINTP_TRI = NINTP_TRI + 1;
         TRINODS(NINTP_TRI) = INOD;
      end

      TRIVERT_IN_FACE(IPT,ITRI) = 1;

   end
   
   % b. Test if Face vertices lay on triangle, including triangle edges:
   for IPF=1:NFVERT
      % Transform back to master Element coordinates
      % location of point i,j in x2-x3 coordinates:
      FD(1:2) = [ FVERT(IAXIS,IPF)-XYEL(IAXIS,NOD3), FVERT(JAXIS,IPF)-XYEL(JAXIS,NOD3) ];
      % Here xi in vec(1) and eta in vec(2)
      VEC(IAXIS) = BODINT_PLANE.AINV(1,1,ITRI)*FD(1) + BODINT_PLANE.AINV(1,2,ITRI)*FD(2);
      VEC(JAXIS) = BODINT_PLANE.AINV(2,1,ITRI)*FD(1) + BODINT_PLANE.AINV(2,2,ITRI)*FD(2);

      % Test for vertex point within triangle, considers Triangle Edges:
      if ( (VEC(IAXIS) >= (0.-GEOMEPS)) && ...
           (VEC(JAXIS) >= (0.-GEOMEPS)) && ...
           (1.-VEC(IAXIS)-VEC(JAXIS) >= (0.-GEOMEPS)) )

         % Insertion add point to intersection list:
         XP(IAXIS:JAXIS) = FVERT(IAXIS:JAXIS,IPF)';
         [NINTP,SIZE_X2X3VERT,X2X3VERT,INOD]=INSERT_POINT_2D(XP,NINTP,SIZE_X2X3VERT,X2X3VERT);

         % Insert sort node to triangles local list
         TRUETHAT = true;
         for INP=1:NINTP_TRI
            if (TRINODS(INP) == INOD)
               TRUETHAT = false;
               break
            end
         end
         if ( TRUETHAT ) % new inod entry on list
            NINTP_TRI = NINTP_TRI + 1;
            TRINODS(NINTP_TRI) = INOD;
         end

         FVERT_IN_TRIANG(IPF,ITRI) = 1;

      end
   end
   
   % Now add face edge - triangle edge intersection points:
   % x2 segments:
  for MYAXIS=IAXIS:JAXIS
      switch(MYAXIS)
         case(IAXIS)
            XIAXIS = IAXIS;
            XJAXIS = JAXIS;
            XIPLNS(LOW_IND:HIGH_IND) = [ X2FMIN, X2FMAX ];
            XJPLNS(LOW_IND:HIGH_IND) = [ X3FMIN, X3FMAX ];
         case(JAXIS)
            XIAXIS = JAXIS;
            XJAXIS = IAXIS;
            XIPLNS(LOW_IND:HIGH_IND) = [ X3FMIN, X3FMAX ];
            XJPLNS(LOW_IND:HIGH_IND) = [ X2FMIN, X2FMAX ];
      end

      for JPL=LOW_IND:HIGH_IND

         XJPLN = XJPLNS(JPL);

         for IPT=1:NTVERT

            XY1(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , TSEGS(NOD1,IPT) )';
            XY2(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , TSEGS(NOD2,IPT) )';

            % Drop if Triangle edge on one side of segment ray:
            MAXXJ = max(XY1(XJAXIS),XY2(XJAXIS));
            MINXJ = min(XY1(XJAXIS),XY2(XJAXIS));
            OUTPLANE1 = ((XJPLN-MAXXJ) > GEOMEPS) || ((MINXJ-XJPLN) > GEOMEPS);
            if ( OUTPLANE1 ); continue; end

            % Also drop if Triangle edge ouside of face edge limits:
            MAXXI = max(XY1(XIAXIS),XY2(XIAXIS));
            MINXI = min(XY1(XIAXIS),XY2(XIAXIS));
            OUTPLANE2 = ((XIPLNS(LOW_IND)-MAXXI) > GEOMEPS) || ((MINXI-XIPLNS(HIGH_IND)) > GEOMEPS);
            if ( OUTPLANE2 ); continue; end

            % Test if segment aligned with xi
            XIALIGNED = ((MAXXJ-MINXJ) < GEOMEPS);
            if ( XIALIGNED ); continue; end % Aligned and on top of xjpln: Intersection points already added.

            % Drop intersections in triangle segment nodes: already added.
            % Compute: dot(plnormal, xyzv - xypl):
            DOT1 = XY1(XJAXIS) - XJPLN;
            DOT2 = XY2(XJAXIS) - XJPLN;

            if ( abs(DOT1) <= GEOMEPS ); continue; end
            if ( abs(DOT2) <= GEOMEPS ); continue; end

            % Finally regular case:
            % Points 1 on one side of x2 segment, point 2 on the other:
            if ( DOT1*DOT2 < 0. )

               % Intersection Point along segment:
               DS    = (XJPLN-XY1(XJAXIS))/(XY2(XJAXIS)-XY1(XJAXIS));
               SVARI = XY1(XIAXIS) + DS*(XY2(XIAXIS)-XY1(XIAXIS));

               OUTSEG= ((XIPLNS(LOW_IND)-SVARI) > -GEOMEPS) || ((SVARI-XIPLNS(HIGH_IND)) > -GEOMEPS);
               if ( OUTSEG ); continue; end

               % Insertion add point to intersection list:
               XP(XIAXIS) = SVARI;
               XP(XJAXIS) = XJPLN;
               [NINTP,SIZE_X2X3VERT,X2X3VERT,INOD]=INSERT_POINT_2D(XP,NINTP,SIZE_X2X3VERT,X2X3VERT);

               % Insert sort node to triangles local list
               TRUETHAT = true;
               for INP=1:NINTP_TRI
                  if (TRINODS(INP) == INOD)
                     TRUETHAT = false;
                     break
                  end
               end
               if (TRUETHAT) % new inod entry on list
                  NINTP_TRI = NINTP_TRI + 1;
                  TRINODS(NINTP_TRI) = INOD;
               end
               continue
            end
         end
      end
  end
  
   if ( NINTP_TRI == 0 ); continue; end

   % Reorder points given normal on x1 direction:
   % Centroid:
   XCEN(IAXIS:JAXIS) = 0.;
   for INTP=1:NINTP_TRI
      XCEN(IAXIS:JAXIS) = XCEN(IAXIS:JAXIS) + X2X3VERT(IAXIS:JAXIS,TRINODS(INTP))';
   end
   XCEN(IAXIS:JAXIS)= XCEN(IAXIS:JAXIS) * NINTP_TRI^(-1.);

   ATANTRI(1:IBM_MAXVERTS_FACE+1) = 1. / GEOMEPS;
   II(1:IBM_MAXVERTS_FACE+1) = IBM_UNDEFINED;
   for INTP=1:NINTP_TRI
      ATTRI = atan2(X2X3VERT(JAXIS,TRINODS(INTP))-XCEN(JAXIS), ...
                    X2X3VERT(IAXIS,TRINODS(INTP))-XCEN(IAXIS)) + pi;
      % Insertion sort:
      for IINS=1:INTP+1
         if (ATTRI < ATANTRI(IINS)); break; end
      end
      % copy from the back:
      for IDUM=INTP+1:-1:IINS+1
         ATANTRI(IDUM) = ATANTRI(IDUM-1);
         II(IDUM)      = II(IDUM-1);
      end
      ATANTRI(IINS) = ATTRI;
      II(IINS)      = INTP;
   end

   % Reorder nodes:
   TRINODS(1:NINTP_TRI) = TRINODS(II(1:NINTP_TRI));

   % Define and Insertion add segments to CFELEM, indseg
   EDGETRI = IBM_UNDEFINED;
   for IEDGE=1:NINTP_TRI-1
        EDGETRI(NOD1:NOD2,IEDGE) = [ TRINODS(IEDGE), TRINODS(IEDGE+1) ]';
   end
   EDGETRI(NOD1:NOD2,NINTP_TRI) = [ TRINODS(NINTP_TRI), TRINODS(1) ]';

   LOCTRI = BODINT_PLANE.INDTRI(1,ITRI);
   LOCBOD = BODINT_PLANE.INDTRI(2,ITRI);

   for IEDGE=1:NINTP_TRI

      if ( EDGETRI(NOD1,IEDGE) == EDGETRI(NOD2,IEDGE) ); continue; end

      % Test if Edge already on list:
      INLIST = false;
      for ISEG=1:NEDGE

         if ( (EDGETRI(NOD1,IEDGE) == CEELEM(NOD1,ISEG)) && ... % same inod1
              (EDGETRI(NOD2,IEDGE) == CEELEM(NOD2,ISEG)) && ... % same inod2
              (LOCBOD              == INDSEG(4,ISEG)) )         % same ibod

            switch(INDSEG(1,ISEG))
               % Only one triangle in list:
               case(1)
                  if ( LOCTRI ~= INDSEG(2,ISEG) )
                     INDSEG(1,ISEG) = 2;
                     INDSEG(3,ISEG) = LOCTRI; % add triangle 2nd.
                  end
                  INLIST = true;
                  break
               % Two triangles in list:
               case(2)
                  if ( (LOCTRI == INDSEG(2,ISEG)) || ...
                       (LOCTRI == INDSEG(3,ISEG)) )
                     INLIST = true;
                     break
                  end
            end
         end
      end

      if ( ~INLIST ) % Edge not in list.
         NEDGE = NEDGE + 1;
         CEELEM(NOD1:NOD2,NEDGE) =  EDGETRI(NOD1:NOD2,IEDGE);

         % Here we have to figure out if segment belongs to a triangles side:
         SEG_IN_SIDE = false;
         for IPT=1:NTVERT

            % Triangle side nodes:
            XY1(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , TSEGS(NOD1,IPT) )';
            XY2(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , TSEGS(NOD2,IPT) )';

            % Segment points:
            XP1(IAXIS:JAXIS) = X2X3VERT(IAXIS:JAXIS,CEELEM(NOD1,NEDGE))';
            XP2(IAXIS:JAXIS) = X2X3VERT(IAXIS:JAXIS,CEELEM(NOD2,NEDGE))';

            VECS(IAXIS:JAXIS)  = XY2(IAXIS:JAXIS) - XY1(IAXIS:JAXIS);
            VECP1(IAXIS:JAXIS) = XP1(IAXIS:JAXIS) - XY1(IAXIS:JAXIS);
            VECP2(IAXIS:JAXIS) = XP2(IAXIS:JAXIS) - XY1(IAXIS:JAXIS);

            CROSSP1 = abs(VECS(IAXIS)*VECP1(JAXIS)-VECS(JAXIS)*VECP1(IAXIS));
            CROSSP2 = abs(VECS(IAXIS)*VECP2(JAXIS)-VECS(JAXIS)*VECP2(IAXIS));

            if ( (CROSSP1+CROSSP2) < GEOMEPS )
               SEG_IN_SIDE = true;
               break
            end
         end
         if ( SEG_IN_SIDE )
            EDGE_TRI = GEOM(LOCBOD).FACE_EDGES(IPT,LOCTRI); % WSTRIED
            VEC3(1)  = GEOM(LOCBOD).EDGE_FACES(1,EDGE_TRI); % WSEDTRI
            VEC3(2)  = GEOM(LOCBOD).EDGE_FACES(2,EDGE_TRI);
            VEC3(3)  = GEOM(LOCBOD).EDGE_FACES(4,EDGE_TRI);
            INDSEG(1:4,NEDGE) = [ VEC3(1), VEC3(2), VEC3(3), LOCBOD ];
         else
            INDSEG(1:4,NEDGE) = [ 1, LOCTRI, 0, LOCBOD ];
         end
      end
   end

end

% Now define cut-edges from solid-solid segments:
for IWSSEG=1:BODINT_PLANE.NSEGS

   NINTP_SEG = 0;
   SEGNODS   = IBM_UNDEFINED;

   SEG(NOD1:NOD2) = BODINT_PLANE.SEGS(NOD1:NOD2,IWSSEG);
   for INOD=NOD1:NOD2
      XYEL(IAXIS:JAXIS,INOD) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] ,SEG(INOD));
   end
   % Cycle if Edges BBOX not intersecting face:
   OUTX2= ((X2FMIN-max(XYEL(IAXIS,NOD1:NOD2))) > GEOMEPS) || ...
          ((min(XYEL(IAXIS,NOD1:NOD2))-X2FMAX) > GEOMEPS); % Segment out of Face in x2 dir
   OUTX3= ((X3FMIN-max(XYEL(JAXIS,NOD1:NOD2))) > GEOMEPS) || ...
          ((min(XYEL(JAXIS,NOD1:NOD2))-X3FMAX) > GEOMEPS); % Segment out of Face in x3 dir
   OUTFACE = OUTX2 || OUTX3;
   if (OUTFACE); continue; end

   % Now define nodes for this CEELEM:
   % a-1. Test if Segments vertices Lay on Faces area, including face boundary:
   for IPT=1:NSVERT
      OUTX2= ((X2FMIN-XYEL(IAXIS,IPT)) > GEOMEPS) || ...
             ((XYEL(IAXIS,IPT)-X2FMAX) > GEOMEPS);  % Triang out of Face in x2 dir
      OUTX3= ((X3FMIN-XYEL(JAXIS,IPT)) > GEOMEPS) || ...
             ((XYEL(JAXIS,IPT)-X3FMAX) > GEOMEPS);  % Triang out of Face in x3 dir
      OUTFACE = OUTX2 || OUTX3;
      if ( OUTFACE ); continue; end

      % Insertion add point to intersection list:
      XP(IAXIS:JAXIS) = XYEL(IAXIS:JAXIS,IPT)';
      [NINTP,SIZE_X2X3VERT,X2X3VERT,INOD]=INSERT_POINT_2D(XP,NINTP,SIZE_X2X3VERT,X2X3VERT);

      % Insert sort node to triangles local list
      TRUETHAT = true;
      for INP=1:NINTP_SEG
         if (SEGNODS(INP) == INOD)
            TRUETHAT = false;
            break
         end
      end
      if ( TRUETHAT ) % new inod entry on list
         NINTP_SEG = NINTP_SEG + 1;
         SEGNODS(NINTP_SEG) = INOD;
      end
   end

   if (NINTP_SEG < 2)
      % b. Now add face edge - SS edge intersection points:
      % x2 segments:
      for MYAXIS=IAXIS:JAXIS
         switch(MYAXIS)
            case(IAXIS)
               XIAXIS = IAXIS;
               XJAXIS = JAXIS;
               XIPLNS(LOW_IND:HIGH_IND) = [ X2FMIN, X2FMAX ];
               XJPLNS(LOW_IND:HIGH_IND) = [ X3FMIN, X3FMAX ];
            case(JAXIS)
               XIAXIS = JAXIS;
               XJAXIS = IAXIS;
               XIPLNS(LOW_IND:HIGH_IND) = [ X3FMIN, X3FMAX ];
               XJPLNS(LOW_IND:HIGH_IND) = [ X2FMIN, X2FMAX ];
         end

         for JPL=LOW_IND:HIGH_IND

            XJPLN = XJPLNS(JPL);

            XY1(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , SEG(NOD1) )';
            XY2(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [ X2AXIS, X3AXIS ] , SEG(NOD2) )';

            % b-1. Drop if Edge on one side of segment ray:
            MAXXJ = max(XY1(XJAXIS),XY2(XJAXIS));
            MINXJ = min(XY1(XJAXIS),XY2(XJAXIS));
            OUTPLANE1 = ((XJPLN-MAXXJ) > GEOMEPS) || ((MINXJ-XJPLN) > GEOMEPS);
            if ( OUTPLANE1 ); continue; end

            % b-2. Also drop if Edge ouside of face edge limits:
            MAXXI = max(XY1(XIAXIS),XY2(XIAXIS));
            MINXI = min(XY1(XIAXIS),XY2(XIAXIS));
            OUTPLANE2 = ((XIPLNS(LOW_IND)-MAXXI) > GEOMEPS) || ((MINXI-XIPLNS(HIGH_IND)) > GEOMEPS);
            if ( OUTPLANE2 ); continue; end

            % Test if segment aligned with xi
            XIALIGNED = ((MAXXJ-MINXJ) < GEOMEPS);
            if ( XIALIGNED ); continue; end % Aligned and on top of xjpln: Intersection points already added.

            % Drop intersections in EDGE nodes: already added.
            % Compute: dot(plnormal, xyzv - xypl):
            DOT1 = XY1(XJAXIS) - XJPLN;
            DOT2 = XY2(XJAXIS) - XJPLN;

            if ( abs(DOT1) <= GEOMEPS ); continue; end
            if ( abs(DOT2) <= GEOMEPS ); continue; end

            % Finally regular case:
            % Points 1 on one side of x2 segment, point 2 on the other:
            if ( DOT1*DOT2 < 0. )

               % Intersection Point along segment:
               DS    = (XJPLN-XY1(XJAXIS))/(XY2(XJAXIS)-XY1(XJAXIS));
               SVARI = XY1(XIAXIS) + DS*(XY2(XIAXIS)-XY1(XIAXIS));

               OUTSEG= ((XIPLNS(LOW_IND)-SVARI) > -GEOMEPS) || ((SVARI-XIPLNS(HIGH_IND)) > -GEOMEPS);
               if ( OUTSEG ); continue; end

               % Insertion add point to intersection list:
               XP(XIAXIS) = SVARI;
               XP(XJAXIS) = XJPLN;
               [NINTP,SIZE_X2X3VERT,X2X3VERT,INOD]=INSERT_POINT_2D(XP,NINTP,SIZE_X2X3VERT,X2X3VERT);

               % Insert sort node to EDGES local list
               TRUETHAT = true;
               for INP=1:NINTP_SEG
                  if (SEGNODS(INP) == INOD)
                     TRUETHAT = false;
                     break
                  end
               end
               if (TRUETHAT) % new inod entry on list
                  NINTP_SEG = NINTP_SEG + 1;
                  SEGNODS(NINTP_SEG) = INOD;
               end
               continue
            end
         end
      end
   end

   if ( (NINTP_SEG < 2) || (SEGNODS(NOD1) == SEGNODS(NOD2)) ); continue; end

   % Test if Edge already on list:
   INLIST = false;
   for ISEG=1:NEDGE

      if ( (SEGNODS(NOD1) == CEELEM(NOD1,ISEG)) && ... % same inod1
           (SEGNODS(NOD2) == CEELEM(NOD2,ISEG)) && ... % same inod2
           (BODINT_PLANE.INDSEG(4,IWSSEG) == INDSEG(4,ISEG)) )     % same ibod

         if (any(BODINT_PLANE.INDSEG(2:3,IWSSEG) == INDSEG(2,ISEG)))
            % Edge already in list, Use SS Edge INDSEG:
            INDSEG(1:4,ISEG) = BODINT_PLANE.INDSEG(1:4,IWSSEG);
            INLIST = true;
            break
         else
            disp('Error in GET_TRIANG_FACE_INT: SS EDGE Triangles not on 2 WS triang list INDSEG.')
         end
      end
   end

   if ( ~INLIST ) % Edge not in list.
      NEDGE = NEDGE + 1;
      CEELEM(NOD1:NOD2,NEDGE) =  SEGNODS(NOD1:NOD2)';
      INDSEG(1:4,NEDGE)       = BODINT_PLANE.INDSEG(1:4,IWSSEG);
   end
end

% Populate XYVERT points array:
XYVERT= X2X3VERT;
NVERT = NINTP;
if (NVERT > 0); INB_FLG = true; end


return