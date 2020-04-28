function [ierr,BODINT_PLANE]=GET_BODINT_PLANE(X1AXIS,X1PLN,INDX1,PLNORMAL,X2AXIS,...
                             X3AXIS,DX2_MIN,DX3_MIN,X2LO,X2HI,X3LO,X3HI,X2FACE,X3FACE,TRI_ONPLANE_ONLY,RAYTRACE_X2_ONLY)
                         
global N_GEOMETRY GEOM IAXIS JAXIS KAXIS IBM_SOLID IBM_GASPHASE
global GEOMEPS LOW_IND HIGH_IND NODS_WSEL MAX_DIM NOD1 NOD2 NOD3 EDG1 EDG2 EDG3
global IBM_MAX_WSTRIANG_SEG DELTA_SEGBIN
global IBM_DELTA_NBCROSS

global X2LO_CELL X2HI_CELL DX2FACE
global X3LO_CELL X3HI_CELL DX3FACE

global FACERT CELLRT
global X2NOC X3NOC

global X3LO_RT X3HI_RT

global BODINT_PLANE3

global XIAXIS XJAXIS XKAXIS


ierr = 1;

% Now allocate BODINT_PLANE:
BODINT_PLANE.X1AXIS  =X1AXIS;
BODINT_PLANE.X2AXIS  =X2AXIS;
BODINT_PLANE.X3AXIS  =X3AXIS;
BODINT_PLANE.X1PLN   =X1PLN;
BODINT_PLANE.GEOMEPS =GEOMEPS;
BODINT_PLANE.NNODS   = 0;
BODINT_PLANE.NSGLS   = 0;
BODINT_PLANE.NSEGS   = 0;
BODINT_PLANE.NTRIS   = 0;
BODINT_PLANE.XYZ     = [];
BODINT_PLANE.NOD_PERM= [];
BODINT_PLANE.SEGS    = [];
BODINT_PLANE.SEGTYPE = [];
BODINT_PLANE.SGLS    = [];
BODINT_PLANE.TRIS    = [];

for IG=1:N_GEOMETRY
    for IWSEL=1:GEOM(IG).N_FACES
       
       % Test low-high vertices of triangle along x1axis vs plane (O(NT) operation):
       if( (GEOM(IG).FACECUBE( LOW_IND,X1AXIS,IWSEL)-X1PLN) > GEOMEPS); continue; end
       if( (X1PLN-GEOM(IG).FACECUBE(HIGH_IND,X1AXIS,IWSEL)) > GEOMEPS); continue; end
               
       if (RAYTRACE_X2_ONLY)
         if( (X3LO_RT-GEOM(IG).FACECUBE(HIGH_IND,X3AXIS,IWSEL)) > GEOMEPS); continue; end
         if( (GEOM(IG).FACECUBE( LOW_IND,X3AXIS,IWSEL)-X3HI_RT) > GEOMEPS); continue; end
       else
         LO_X2_TEST = (X2FACE(X2LO)-GEOM(IG).FACECUBE(HIGH_IND,X2AXIS,IWSEL)) > GEOMEPS;
         LO_X3_TEST = (X3FACE(X3LO)-GEOM(IG).FACECUBE(HIGH_IND,X3AXIS,IWSEL)) > GEOMEPS;
         if(  LO_X2_TEST && LO_X3_TEST ); continue; end
         HI_X2_TEST = (GEOM(IG).FACECUBE( LOW_IND,X2AXIS,IWSEL)-X2FACE(X2HI)) > GEOMEPS;
         if(  HI_X2_TEST && LO_X3_TEST ); continue; end
         HI_X3_TEST = (GEOM(IG).FACECUBE( LOW_IND,X3AXIS,IWSEL)-X3FACE(X3HI)) > GEOMEPS;
         if(  LO_X2_TEST && HI_X3_TEST ); continue; end
         if(  HI_X2_TEST && HI_X3_TEST ); continue; end
       end
       
       WSELEM(NOD1:NOD3) = GEOM(IG).FACES(NODS_WSEL*(IWSEL-1)+1:NODS_WSEL*IWSEL);
       % Triangles NODES coordinates:
       for INOD=NOD1:NOD3
          XYZV(IAXIS:KAXIS,INOD) = GEOM(IG).VERTS(MAX_DIM*(WSELEM(INOD)-1)+1:MAX_DIM*WSELEM(INOD));
       end

      % Compute simplified dot(PLNORMAL,XYZV-XYZPLANE):
      DOT1 = XYZV(X1AXIS,NOD1) - X1PLN;
      DOT2 = XYZV(X1AXIS,NOD2) - X1PLN;
      DOT3 = XYZV(X1AXIS,NOD3) - X1PLN;
      if ( abs(DOT1) <= GEOMEPS ); DOT1 = 0.; end
      if ( abs(DOT2) <= GEOMEPS ); DOT2 = 0.; end
      if ( abs(DOT3) <= GEOMEPS ); DOT3 = 0.; end

      
      % Test if IWSEL lays in X1PLN:
      if ( (abs(DOT1)+abs(DOT2)+abs(DOT3)) == 0. )

         % Force nodes location in X1PLN plane:
         XYZV(X1AXIS,NOD1:NOD3) = X1PLN;

         % Index to point 1 of triangle in BODINT_PLANE.XYZ list:
         [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZV(IAXIS:KAXIS,NOD1));

         % Index to point 2 of triangle in BODINT_PLANE.XYZ list:
         [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZV(IAXIS:KAXIS,NOD2));

         % Index to point 3 of triangle in BODINT_PLANE.XYZ list:
         [IND_P(NOD3),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZV(IAXIS:KAXIS,NOD3));

         % Do we need to test if we already have this triangle on
         % the list? Shouldn't unless repeated -> Possibility for
         % zero thickness.
         NTRIS = BODINT_PLANE.NTRIS + 1;
         BODINT_PLANE.NTRIS = NTRIS;
         BODINT_PLANE.TRIS(NOD1:NOD3,NTRIS) = IND_P';
         BODINT_PLANE.INDTRI(1:2,NTRIS) = [ IWSEL, IG ]';

         continue % Next WSELEM

      end

      % Test if we are looking for intersection triangles only:
      if (~TRI_ONPLANE_ONLY)
         % Case a: Typical intersections:
         % Points 1,2 on on side of plane, point 3 on the other:
         if ( ((DOT1 > 0.) && (DOT2 > 0.) && (DOT3 < 0.)) || ...
              ((DOT1 < 0.) && (DOT2 < 0.) && (DOT3 > 0.)) )

            % Line 1, from node 2 to 3:
            LN1(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD2);
            LN1(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

            [XYZ_INT1,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN1);

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Line 2, from node 1 to 3:
            LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
            LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

            [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

            % Index to XYZ_INT2:
            [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

            % Now add segment:
            NSEGS = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = NSEGS;
            if ( DOT1 > 0. ) % First case, counterclockwise p1 to p2
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
            else
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
            end
            BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
            BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

            continue % Next WSELEM

         end
         % Points 2,3 on one side of plane, point 1 on the other:
         if ( ((DOT2 > 0.) && (DOT3 > 0.) && (DOT1 < 0.)) || ...
              ((DOT2 < 0.) && (DOT3 < 0.) && (DOT1 > 0.)) )

              % Line 1, from node 1 to 2:
              LN1(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
              LN1(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD2);

              [XYZ_INT1,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN1);

              % Index to XYZ_INT1:
              [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

              % Line 2, from node 1 to 3:
              LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
              LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

              [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

              % Index to XYZ_INT2:
              [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

              % Now add segment:
              NSEGS = BODINT_PLANE.NSEGS + 1;
              BODINT_PLANE.NSEGS = NSEGS;
              if ( DOT2 > 0. ) % Second case, counterclockwise p2 to p1
                 BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
              else
                 BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
              end
              BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
              BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

              continue % Next WSELEM

         end
         % Points 1,3 on one side of plane, point 2 on the other:
         if ( ((DOT1 > 0.) && (DOT3 > 0.) && (DOT2 < 0.)) || ...
              ((DOT1 < 0.) && (DOT3 < 0.) && (DOT2 > 0.)) )

              % Line 1, from node 1 to 2:
              LN1(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
              LN1(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD2);

              [XYZ_INT1,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN1);

              % Index to XYZ_INT1:
              [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

              % Line 2, from node 2 to 3:
              LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD2);
              LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

              [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

              % Index to XYZ_INT2:
              [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

              % Now add segment:
              NSEGS = BODINT_PLANE.NSEGS + 1;
              BODINT_PLANE.NSEGS = NSEGS;
              if ( DOT1 > 0. ) % Third case, counterclockwise p1 to p2
                 BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
              else
                 BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
              end
              BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
              BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

              continue % Next WSELEM

         end

         % Case b: only one point intersection. They will be used to define
         % Solid vertex points in case of coincidence.
         % Point 1 is on the plane:
         if ( (DOT1 == 0.) && ( ((DOT2 > 0.) && (DOT3 > 0.)) || ...
                                ((DOT2 < 0.) && (DOT3 < 0.)) ) )

            % First node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD1); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Add index to singles:
            % Find if oriented segment is in list:
            INLIST = false;
            for ISGL=1:BODINT_PLANE.NSGLS
               if (BODINT_PLANE.SGLS(NOD1,ISGL) == IND_P(NOD1))
                  INLIST = true;
                  break
               end
            end
            if (~INLIST)
               ISGL = BODINT_PLANE.NSGLS + 1;
               BODINT_PLANE.NSGLS = ISGL;
               BODINT_PLANE.SGLS(NOD1,ISGL) = IND_P(NOD1);
            end

            continue % Next WSELEM

         end
         % Point 2 is on the plane:
         if ( (DOT2 == 0.) && ( ((DOT1 > 0.) && (DOT3 > 0.)) || ...
                                ((DOT1 < 0.) && (DOT3 < 0.)) ) )

            % Second node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD2); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Add index to singles:
            % Find if oriented segment is in list:
            INLIST = false;
            for ISGL=1:BODINT_PLANE.NSGLS
               if (BODINT_PLANE.SGLS(NOD1,ISGL) == IND_P(NOD1))
                  INLIST = true;
                  break
               end
            end
            if (~INLIST)
               ISGL = BODINT_PLANE.NSGLS + 1;
               BODINT_PLANE.NSGLS = ISGL;
               BODINT_PLANE.SGLS(NOD1,ISGL) = IND_P(NOD1);
            end

            continue % Next WSELEM

         end
         % Point 3 is on the plane:
         if ( (DOT3 == 0.) && ( ((DOT1 > 0.) && (DOT2 > 0.)) || ...
                                ((DOT1 < 0.) && (DOT2 < 0.)) ) )

            % Third node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD3); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Add index to singles:
            % Find if single element is in list:
            INLIST = false;
            for ISGL=1:BODINT_PLANE.NSGLS
               if (BODINT_PLANE.SGLS(NOD1,ISGL) == IND_P(NOD1))
                  INLIST = true;
                  break
               end
            end
            if (~INLIST)
               ISGL = BODINT_PLANE.NSGLS + 1;
               BODINT_PLANE.NSGLS = ISGL;
               BODINT_PLANE.SGLS(NOD1,ISGL) = IND_P(NOD1);
            end

            continue % Next WSELEM

         end

         % Case c: one node is part of the intersection:
         % Node 1 is in the plane:
         if ( (DOT1 == 0.) && ( ((DOT2 > 0.) && (DOT3 < 0.)) || ...
                                ((DOT2 < 0.) && (DOT3 > 0.)) ) )

            % First node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD1); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Line 2, from node 2 to 3:
            LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD2);
            LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

            [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

            % Index to XYZ_INT2:
            [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

            % Now add segment:
            NSEGS = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = NSEGS;
            if ( DOT2 > 0. ) % Second case, counterclockwise p2 to p1
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
            else
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
            end
            BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
            BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

            continue % Next WSELEM

         end
         % Node 2 is in the plane:
         if ( (DOT2 == 0.) && ( ((DOT1 > 0.) && (DOT3 < 0.)) || ...
                                ((DOT1 < 0.) && (DOT3 > 0.)) ) )

            % Second node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD2); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Line 2, from node 1 to 3:
            LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
            LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD3);

            [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

            % Index to XYZ_INT2:
            [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

            % Now add segment:
            NSEGS = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = NSEGS;
            if ( DOT1 > 0. )
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
            else
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
            end
            BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
            BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

            continue % Next WSELEM

         end
         % Node 3 is in the plane:
         if ( (DOT3 == 0.) && ( ((DOT1 > 0.) && (DOT2 < 0.)) || ...
                                ((DOT1 < 0.) && (DOT2 > 0.)) ) )

            % Third node is an intersection point:
            XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD3); XYZ_INT1(X1AXIS) = X1PLN;

            % Index to XYZ_INT1:
            [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

            % Line 2, from node 1 to 2:
            LN2(IAXIS:KAXIS,NOD1) = XYZV(IAXIS:KAXIS,NOD1);
            LN2(IAXIS:KAXIS,NOD2) = XYZV(IAXIS:KAXIS,NOD2);

            [XYZ_INT2,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LN2);

            % Index to XYZ_INT2:
            [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

            % Now add segment:
            NSEGS = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = NSEGS;
            if ( DOT1 > 0. )
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD2), IND_P(NOD1) ]';
            else
               BODINT_PLANE.SEGS(NOD1:NOD2,NSEGS) = [ IND_P(NOD1), IND_P(NOD2) ]';
            end
            BODINT_PLANE.INDSEG(1:4,NSEGS) = [ 1, IWSEL, 0, IG ]';
            BODINT_PLANE.SEGTYPE(1:2,NSEGS)= [ IBM_SOLID, IBM_GASPHASE ]';

            continue % Next WSELEM

         end
      end

      % Case D: A triangle segment is in the plane.
      % Intersection is line 1-2:
      if ( (DOT1 == 0.) && (DOT2 == 0.) )

         % First node:
         XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD1); XYZ_INT1(X1AXIS) = X1PLN;

         % Index to XYZ_INT1:
         [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

         % Second node:
         XYZ_INT2 = XYZV(IAXIS:KAXIS,NOD2); XYZ_INT2(X1AXIS) = X1PLN;

         % Index to XYZ_INT2:
         [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

         % Set oriented segment regarding plane:
         if ( DOT3 > 0. )
            SEG(NOD1:NOD2) = [ IND_P(NOD1), IND_P(NOD2) ];
         else
            SEG(NOD1:NOD2) = [ IND_P(NOD2), IND_P(NOD1) ];
         end
         % Find if oriented segment is in list:
         EDGE_TRI = GEOM(IG).FACE_EDGES(EDG1,IWSEL); % 1st edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
         VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
         VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
         VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
         INLIST = false;
         for ISEG=1:BODINT_PLANE.NSEGS
            SEG_FLG = ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD1))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD2))) || ...
                      ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD2))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD1)));
            if (  SEG_FLG && ...
                 (BODINT_PLANE.INDSEG(2,ISEG)  == VEC3(2))   && ...
                 (BODINT_PLANE.INDSEG(3,ISEG)  == VEC3(3))   && ...
                 (BODINT_PLANE.INDSEG(4,ISEG)  == IG))
               INLIST = true;
               break
            end
         end
         if (~INLIST)
            ISEG = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = ISEG;
            BODINT_PLANE.SEGS(NOD1:NOD2,ISEG) = SEG;
%             EDGE_TRI = GEOM(IG).FACE_EDGES(EDG1,IWSEL); % 1st edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
%             VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
%             VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
%             VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
            BODINT_PLANE.INDSEG(1:4,ISEG) = [ VEC3(1), VEC3(2), VEC3(3), IG ]';
         end

         continue % Next WSELEM

      end
      % Intersection is line 2-3:
      if ( (DOT2 == 0.) && (DOT3 == 0.) )

         % Second node:
         XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD2); XYZ_INT1(X1AXIS) = X1PLN;

         % Index to XYZ_INT1:
         [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

         % Third node:
         XYZ_INT2 = XYZV(IAXIS:KAXIS,NOD3); XYZ_INT2(X1AXIS) = X1PLN;

         % Index to XYZ_INT2:
         [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

         % Set oriented segment regarding plane:
         if ( DOT1 > 0. )
            SEG(NOD1:NOD2) = [ IND_P(NOD1), IND_P(NOD2) ];
         else
            SEG(NOD1:NOD2) = [ IND_P(NOD2), IND_P(NOD1) ];
         end
         % Find if oriented segment is in list:
         EDGE_TRI = GEOM(IG).FACE_EDGES(EDG2,IWSEL); % 2nd edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
         VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
         VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
         VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
         INLIST = false;
         for ISEG=1:BODINT_PLANE.NSEGS
            SEG_FLG = ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD1))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD2))) || ...
                      ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD2))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD1)));
            if (  SEG_FLG && ...
                 (BODINT_PLANE.INDSEG(2,ISEG)  == VEC3(2))   && ...
                 (BODINT_PLANE.INDSEG(3,ISEG)  == VEC3(3))   && ...
                 (BODINT_PLANE.INDSEG(4,ISEG)  == IG))
               INLIST = true;
               break
            end
         end
         if (~INLIST)
            ISEG = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = ISEG;
            BODINT_PLANE.SEGS(NOD1:NOD2,ISEG) = SEG;
%             EDGE_TRI = GEOM(IG).FACE_EDGES(EDG2,IWSEL); % 2nd edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
%             VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
%             VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
%             VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
            BODINT_PLANE.INDSEG(1:4,ISEG) = [ VEC3(1), VEC3(2), VEC3(3), IG ]';
         end

         continue % Next WSELEM

      end
      % Intersection is line 3-1:
      if ( (DOT3 == 0.) && (DOT1 == 0.) )

         % Third node:
         XYZ_INT1 = XYZV(IAXIS:KAXIS,NOD3); XYZ_INT1(X1AXIS) = X1PLN;

         % Index to XYZ_INT1:
         [IND_P(NOD1),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT1);

         % First node:
         XYZ_INT2 = XYZV(IAXIS:KAXIS,NOD1); XYZ_INT2(X1AXIS) = X1PLN;

         % Index to XYZ_INT2:
         [IND_P(NOD2),BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ_INT2);

         % Set oriented segment regarding plane:
         if ( DOT2 > 0. )
            SEG(NOD1:NOD2) = [ IND_P(NOD1), IND_P(NOD2) ];
         else
            SEG(NOD1:NOD2) = [ IND_P(NOD2), IND_P(NOD1) ];
         end
         % Find if oriented segment is in list:
         EDGE_TRI = GEOM(IG).FACE_EDGES(EDG3,IWSEL); % 3rd edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
         VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
         VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
         VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
         INLIST = false;
         for ISEG=1:BODINT_PLANE.NSEGS
            SEG_FLG = ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD1))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD2))) || ...
                      ((BODINT_PLANE.SEGS(NOD1,ISEG) == SEG(NOD2))  && ...
                       (BODINT_PLANE.SEGS(NOD2,ISEG) == SEG(NOD1)));
            if (  SEG_FLG && ...
                 (BODINT_PLANE.INDSEG(2,ISEG)  == VEC3(2))   && ...
                 (BODINT_PLANE.INDSEG(3,ISEG)  == VEC3(3))   && ...
                 (BODINT_PLANE.INDSEG(4,ISEG)  == IG))
               INLIST = true;
               break
            end
         end
         if (~INLIST)
            ISEG = BODINT_PLANE.NSEGS + 1;
            BODINT_PLANE.NSEGS = ISEG;
            BODINT_PLANE.SEGS(NOD1:NOD2,ISEG) = SEG;
%             EDGE_TRI = GEOM(IG).FACE_EDGES(EDG3,IWSEL); % 3rd edge: Ed1 NOD1-NOD2, Ed2 NOD2-NOD3, Ed3 NOD3-NOD1.
%             VEC3(1) = GEOM(IG).EDGE_FACES(1,EDGE_TRI);
%             VEC3(2) = GEOM(IG).EDGE_FACES(2,EDGE_TRI);
%             VEC3(3) = GEOM(IG).EDGE_FACES(4,EDGE_TRI);
            BODINT_PLANE.INDSEG(1:4,ISEG) = [ VEC3(1), VEC3(2), VEC3(3), IG ]';
         end

         continue % Next WSELEM

      end

      % If you get to this point -> you have a problem:
      if (~TRI_ONPLANE_ONLY) 
          disp(['Error GET_BODINT_PLANE: Missed wet surface Triangle =' num2str(IWSEL)])
      end   
       
    end
end    


% Next step is to Test triangles sides normals on plane against the obtained
% segments normals. If two identical segments found contain oposite
% normals, drop the segment in BODINT_PLANE.SEGS:
if ( BODINT_PLANE.NTRIS > 0 )

   for ITRI=1:BODINT_PLANE.NTRIS

      % Triang conectivities:
      ELEM(NOD1:NOD3) = BODINT_PLANE.TRIS(NOD1:NOD3,ITRI);

      % Coordinates in x2, x3 directions:
      X2X3(IAXIS,NOD1:NOD3) = [ BODINT_PLANE.XYZ(X2AXIS,ELEM(NOD1)), ...
                                BODINT_PLANE.XYZ(X2AXIS,ELEM(NOD2)), ...
                                BODINT_PLANE.XYZ(X2AXIS,ELEM(NOD3)) ];
      X2X3(JAXIS,NOD1:NOD3) = [ BODINT_PLANE.XYZ(X3AXIS,ELEM(NOD1)), ...
                                BODINT_PLANE.XYZ(X3AXIS,ELEM(NOD2)), ...
                                BODINT_PLANE.XYZ(X3AXIS,ELEM(NOD3)) ];

      % Test Area sign, if -ve switch node order:
      AREALOC = 0.5*(X2X3(IAXIS,NOD1)*X2X3(JAXIS,NOD2) - X2X3(IAXIS,NOD2)*X2X3(JAXIS,NOD1) + ...
                     X2X3(IAXIS,NOD2)*X2X3(JAXIS,NOD3) - X2X3(IAXIS,NOD3)*X2X3(JAXIS,NOD2) + ...
                     X2X3(IAXIS,NOD3)*X2X3(JAXIS,NOD1) - X2X3(IAXIS,NOD1)*X2X3(JAXIS,NOD3));
      if (AREALOC < 0.)
         ISEG    = ELEM(3);
         ELEM(3) = ELEM(2);
         ELEM(2)  =   ISEG;
      end

      % Now corresponding segments, ordered normal outside of plane x2-x3.
      EDGES(NOD1:NOD2,1) = [ ELEM(1), ELEM(2) ]'; % edge 1.
      EDGES(NOD1:NOD2,2) = [ ELEM(2), ELEM(3) ]'; % edge 2.
      EDGES(NOD1:NOD2,3) = [ ELEM(3), ELEM(1) ]';

      % Now Test against segments, Beast approach:
      for IEDGE=1:3
         for ISEG=1:BODINT_PLANE.NSEGS
            if ( (BODINT_PLANE.SEGS(NOD1,ISEG) == EDGES(NOD2,IEDGE)) && ...
                 (BODINT_PLANE.SEGS(NOD2,ISEG) == EDGES(NOD1,IEDGE)) ) % Edge normals
                                                                              % oriented in opposite dirs.
               % Set to SOLID SOLID segtype from BODINT_PLANE.SEGS
               BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG)=[ IBM_SOLID, IBM_SOLID ]';

            end
         end
      end

   end
end

% For segments that are related to 2 Wet Surface triangles, test if they are of type GG or SS:
for ISEG=1:BODINT_PLANE.NSEGS
    if (BODINT_PLANE.INDSEG(1,ISEG) > 1) % Related to 2 WS triangles:

       SEG(NOD1:NOD2) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);

       % Segment nodes positions:
       XP1(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [X2AXIS,X3AXIS] ,SEG(NOD1));
       XP2(IAXIS:JAXIS) = BODINT_PLANE.XYZ( [X2AXIS,X3AXIS] ,SEG(NOD2));

       % Unit normal versor along x2p (axis directed from NOD2 to NOD1):
       NMTX2P = sqrt( (XP1(IAXIS)-XP2(IAXIS))^2. + (XP1(JAXIS)-XP2(JAXIS))^2. );
       TX2P(IAXIS:JAXIS) = (XP1(IAXIS:JAXIS)-XP2(IAXIS:JAXIS)) * NMTX2P^(-1.);
       % Versor along x3p.
       TX3P(IAXIS:JAXIS) = [ -TX2P(JAXIS), TX2P(IAXIS) ];

       % Now related WS triangles centroids:
       IWSEL1 = BODINT_PLANE.INDSEG(2,ISEG);
       IWSEL2 = BODINT_PLANE.INDSEG(3,ISEG);
       IG     = BODINT_PLANE.INDSEG(4,ISEG);

       % Centroid of WS elem 1:
       ELEM1(NOD1:NOD3)  =   GEOM(IG).FACES(NODS_WSEL*(IWSEL1-1)+1:NODS_WSEL*IWSEL1);
       XYZ1(IAXIS:KAXIS) = ( GEOM(IG).VERTS(MAX_DIM*(ELEM1(NOD1)-1)+1:MAX_DIM*ELEM1(NOD1)) + ...
                             GEOM(IG).VERTS(MAX_DIM*(ELEM1(NOD2)-1)+1:MAX_DIM*ELEM1(NOD2)) + ...
                             GEOM(IG).VERTS(MAX_DIM*(ELEM1(NOD3)-1)+1:MAX_DIM*ELEM1(NOD3)) ) / 3.;
       NXYZ1(IAXIS:KAXIS)= GEOM(IG).FACES_NORMAL(IAXIS:KAXIS,IWSEL1);
       % Normal versor in x3p-x1 direction:
       NX3P1 = TX3P(IAXIS)*NXYZ1(X2AXIS) + TX3P(JAXIS)*NXYZ1(X3AXIS);
       N1(IAXIS:JAXIS) = [ NX3P1, NXYZ1(X1AXIS) ];
       NMNL = sqrt( N1(IAXIS)^2. + N1(JAXIS)^2. );
       N1 = N1 * NMNL^(-1.);

       % Centroid of WS elem 2:
       ELEM2(NOD1:NOD3)  = GEOM(IG).FACES(NODS_WSEL*(IWSEL2-1)+1:NODS_WSEL*IWSEL2);
       XYZ2(IAXIS:KAXIS) = ( GEOM(IG).VERTS(MAX_DIM*(ELEM2(NOD1)-1)+1:MAX_DIM*ELEM2(NOD1)) + ...
                             GEOM(IG).VERTS(MAX_DIM*(ELEM2(NOD2)-1)+1:MAX_DIM*ELEM2(NOD2)) + ...
                             GEOM(IG).VERTS(MAX_DIM*(ELEM2(NOD3)-1)+1:MAX_DIM*ELEM2(NOD3)) ) / 3.;
       NXYZ2(IAXIS:KAXIS)= GEOM(IG).FACES_NORMAL(IAXIS:KAXIS,IWSEL2);
       % Normal versor in x3p-x1 direction:
       NX3P2 = TX3P(IAXIS)*NXYZ2(X2AXIS) + TX3P(JAXIS)*NXYZ2(X3AXIS);
       N2(IAXIS:JAXIS) = [ NX3P2, NXYZ2(X1AXIS) ];
       NMNL = sqrt( N2(IAXIS)^2. + N2(JAXIS)^2. );
       N2 = N2 * NMNL^(-1.);

       % Define points in plane x3p-x1:
       % vertex point:
       X3PVERT = TX3P(IAXIS)*XP1(IAXIS) + TX3P(JAXIS)*XP1(JAXIS);
       PVERT(IAXIS:JAXIS) = [ X3PVERT, X1PLN ];
       % First triangle centroid:
       X3P1 = TX3P(IAXIS)*XYZ1(X2AXIS) + TX3P(JAXIS)*XYZ1(X3AXIS);
       P1CEN(IAXIS:JAXIS) = [ X3P1, XYZ1(X1AXIS) ];
       % Second triangle centroid:
       X3P2 = TX3P(IAXIS)*XYZ2(X2AXIS) + TX3P(JAXIS)*XYZ2(X3AXIS);
       P2CEN(IAXIS:JAXIS) = [ X3P2, XYZ2(X1AXIS) ];

       VCT(1:2) = 0;
       PCT(IAXIS:JAXIS,1:2) = 0.;

       % Segment on triangle 1:
       V1(IAXIS:JAXIS) = P1CEN(IAXIS:JAXIS) - PVERT(IAXIS:JAXIS);
       CRSSNV = N1(IAXIS)*V1(JAXIS) - N1(JAXIS)*V1(IAXIS);
       if (CRSSNV > 0.)
           % v1 stays as is, and is second segment:
           VEC(IAXIS:JAXIS,2) = V1(IAXIS:JAXIS)';
           PCT(IAXIS:JAXIS,2) = P1CEN(IAXIS:JAXIS)';
           VCT(2) = 1;
       else
           % -v1 is the first segment:
           VEC(IAXIS:JAXIS,1) = -V1(IAXIS:JAXIS)';
           PCT(IAXIS:JAXIS,1) = P1CEN(IAXIS:JAXIS)';
           VCT(1) = 1;
       end

       % Segment on triangle 2:
       V2(IAXIS:JAXIS) = P2CEN(IAXIS:JAXIS) - PVERT(IAXIS:JAXIS);
       CRSSNV = N2(IAXIS)*V2(JAXIS) - N2(JAXIS)*V2(IAXIS);
       if (CRSSNV > 0.)
           % v2 stays as is, and is second segment:
           VEC(IAXIS:JAXIS,2) = V2(IAXIS:JAXIS)';
           PCT(IAXIS:JAXIS,2) = P2CEN(IAXIS:JAXIS)';
           VCT(2) = 1;
       else
           % -v2 is the first segment:
           VEC(IAXIS:JAXIS,1) = -V2(IAXIS:JAXIS)';
           PCT(IAXIS:JAXIS,1) = P2CEN(IAXIS:JAXIS)';
           VCT(1) = 1;
       end

       if ( (VCT(1) == 0) || (VCT(2) == 0) )
          disp(['Error GET_BODINT_PLANE: One component of vct == 0.'])
       end

       % Cross product of v1 and v2 gives magnitude along x2p axis:
       CTST = VEC(IAXIS,1)*VEC(JAXIS,2) - VEC(JAXIS,1)*VEC(IAXIS,2);

       % Now tests:
       % Start with SOLID GASPHASE  definition for segtype:
       BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG) = [ IBM_SOLID, IBM_GASPHASE ]';

       % Test for SOLID SOLID condition:
       if ( ((PCT(JAXIS,1)-X1PLN) > -GEOMEPS) &&  ...
            ((PCT(JAXIS,2)-X1PLN) > -GEOMEPS) && (CTST < GEOMEPS) )
           BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG) = [ IBM_SOLID, IBM_SOLID ]';
           continue
       elseif (((PCT(JAXIS,1)-X1PLN) < GEOMEPS) && ...
               ((PCT(JAXIS,2)-X1PLN) < GEOMEPS) && (CTST < GEOMEPS) )
           BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG) = [ IBM_SOLID, IBM_SOLID ]';
           continue
       end

       % Test for GASPHASE GASPHASE condition:
       if ( ((PCT(JAXIS,1)-X1PLN) > GEOMEPS) &&  ...
            ((PCT(JAXIS,2)-X1PLN) > GEOMEPS) && (CTST > GEOMEPS) )
            BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG) = [ IBM_GASPHASE, IBM_GASPHASE ]';
            continue
       elseif (((PCT(JAXIS,1)-X1PLN) < -GEOMEPS) &&  ...
               ((PCT(JAXIS,2)-X1PLN) < -GEOMEPS) && (CTST > GEOMEPS) )
            BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG) = [ IBM_GASPHASE, IBM_GASPHASE ]';
            continue
       end

    end
end


% For the time being, as BODINT_PLANE is used to create Cartesian face cut-faces
% We eliminate from the list the SEGTYPE=[SOLID SOLID] segments:
SEGAUX    = zeros(NOD2,BODINT_PLANE.NSEGS);
INDSEGAUX = zeros(IBM_MAX_WSTRIANG_SEG+2,BODINT_PLANE.NSEGS);
SEGTYPEAUX= zeros(NOD2,BODINT_PLANE.NSEGS);

ISEG_NEW = 0;
if (~TRI_ONPLANE_ONLY)
   for ISEG=1:BODINT_PLANE.NSEGS
       if ( (BODINT_PLANE.SEGTYPE(NOD1,ISEG) == IBM_SOLID) && ...
            (BODINT_PLANE.SEGTYPE(NOD2,ISEG) == IBM_SOLID) ); continue; end

          ISEG_NEW = ISEG_NEW + 1;
          SEGAUX(NOD1:NOD2,ISEG_NEW) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
          INDSEGAUX(1:IBM_MAX_WSTRIANG_SEG+2,ISEG_NEW) = ...
             BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,ISEG);
          SEGTYPEAUX(NOD1:NOD2,ISEG_NEW) = BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG);
   end
else
   for ISEG=1:BODINT_PLANE.NSEGS
       if ( (BODINT_PLANE.SEGTYPE(NOD1,ISEG) == IBM_SOLID) && ...
            (BODINT_PLANE.SEGTYPE(NOD2,ISEG) == IBM_SOLID) )

          ISEG_NEW = ISEG_NEW + 1;
          SEGAUX(NOD1:NOD2,ISEG_NEW) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
          INDSEGAUX(1:IBM_MAX_WSTRIANG_SEG+2,ISEG_NEW) = ...
             BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,ISEG);
          SEGTYPEAUX(NOD1:NOD2,ISEG_NEW) = BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG);
       end
   end
end

BODINT_PLANE.NSEGS = ISEG_NEW;
BODINT_PLANE.SEGS(NOD1:NOD2,1:ISEG_NEW) = SEGAUX(NOD1:NOD2,1:ISEG_NEW);
BODINT_PLANE.INDSEG(1:IBM_MAX_WSTRIANG_SEG+2,1:ISEG_NEW) = INDSEGAUX(1:IBM_MAX_WSTRIANG_SEG+2,1:ISEG_NEW);
BODINT_PLANE.SEGTYPE(NOD1:NOD2,1:ISEG_NEW) = SEGTYPEAUX(NOD1:NOD2,1:ISEG_NEW);

if (BODINT_PLANE.NSEGS == 0); return; end
if (TRI_ONPLANE_ONLY); return; end

% Segments Crossings fields:
BODINT_PLANE.NBCROSS = zeros(1,BODINT_PLANE.NSEGS);
BODINT_PLANE.SVAR    = -ones(IBM_DELTA_NBCROSS,BODINT_PLANE.NSEGS);
BODINT_PLANE.BOX = zeros(HIGH_IND, MAX_DIM);
BODINT_PLANE.BOX(LOW_IND, X2AXIS) = min(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))-10.*GEOMEPS;
BODINT_PLANE.BOX(HIGH_IND,X2AXIS) = max(BODINT_PLANE.XYZ(X2AXIS,1:BODINT_PLANE.NNODS))+10.*GEOMEPS;
BODINT_PLANE.BOX(LOW_IND, X3AXIS) = min(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))-10.*GEOMEPS;
BODINT_PLANE.BOX(HIGH_IND,X3AXIS) = max(BODINT_PLANE.XYZ(X3AXIS,1:BODINT_PLANE.NNODS))+10.*GEOMEPS;
if (RAYTRACE_X2_ONLY)
   AXIS = X3AXIS;
   BODINT_PLANE.TBAXIS(AXIS).DELBIN = BODINT_PLANE.BOX(HIGH_IND,AXIS)-BODINT_PLANE.BOX(LOW_IND,AXIS);
   BODINT_PLANE.TBAXIS(AXIS).N_BINS = 1;
   IBIN=1;
   BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).X_LOW  = BODINT_PLANE.BOX( LOW_IND,AXIS);
   BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).X_HIGH = BODINT_PLANE.BOX(HIGH_IND,AXIS);
   BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).NTL    = BODINT_PLANE.NSEGS;
   BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).TRI_LIST=[1:BODINT_PLANE.NSEGS];
   return; 
end

% Initialize nbcross with segment nodes locations:
% Add segment ends as crossings:
SEGS_NODE=zeros(1,BODINT_PLANE.NNODS);
MEAN_SLEN=0.;
for ISEG=1:BODINT_PLANE.NSEGS

   % End nodes to cross:
   SEG(NOD1:NOD2) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
   
   if(any(BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG)~=IBM_GASPHASE))
      SEGS_NODE(SEG(NOD1)) = SEGS_NODE(SEG(NOD1))+1;
      SEGS_NODE(SEG(NOD2)) = SEGS_NODE(SEG(NOD2))+1;
   end
   
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD1));
   XYZ2(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD2));

   % Is segment aligned with x3 direction?
   BODINT_PLANE.X3ALIGNED(ISEG) = (abs(XYZ2(X2AXIS)-XYZ1(X2AXIS)) < GEOMEPS);
   % Is segment aligned with x2 rays?:
   BODINT_PLANE.X2ALIGNED(ISEG) = (abs(XYZ2(X3AXIS)-XYZ1(X3AXIS)) < GEOMEPS);   
   
   % x2_x3 of segment point 1:
   X2_1 = XYZ1(X2AXIS); X3_1 = XYZ1(X3AXIS);
   % x2_x3 of segment point 2:
   X2_2 = XYZ2(X2AXIS); X3_2 = XYZ2(X3AXIS);

   % Segment length:
   SLEN = sqrt( (X2_2-X2_1)^2. + (X3_2-X3_1)^2. );
   MEAN_SLEN = MEAN_SLEN + SLEN;

   % First node:
   SBOD = 0.;
   % Add crossing to BODINT_PLANE:
   NBCROSS = BODINT_PLANE.NBCROSS(ISEG) + 1;
   BODINT_PLANE.NBCROSS(ISEG) = NBCROSS;
   BODINT_PLANE.SVAR(NBCROSS,ISEG) = SBOD;

   % Second node:
   SBOD = SLEN;
   % Add crossing to BODINT_PLANE:
   NBCROSS = BODINT_PLANE.NBCROSS(ISEG) + 1;
   BODINT_PLANE.NBCROSS(ISEG) = NBCROSS;
   BODINT_PLANE.SVAR(NBCROSS,ISEG) = SBOD;

end

% Spread Segments in BINs in the x2-x3 directions:
MEAN_SLEN = MEAN_SLEN / BODINT_PLANE.NSEGS;
VAXIS = [ X2AXIS, X3AXIS ];
for I = 1:2
   AXIS = VAXIS(I);
   LXI  = BODINT_PLANE.BOX(HIGH_IND,AXIS)-BODINT_PLANE.BOX(LOW_IND,AXIS);
   if (BODINT_PLANE.NSEGS < 100)
      BODINT_PLANE.TBAXIS(AXIS).N_BINS = max( 1,ceil(LXI/(MEAN_SLEN)));
   else
      BODINT_PLANE.TBAXIS(AXIS).N_BINS = max(10,ceil(LXI/(MEAN_SLEN)));
   end
   
   % Allocate TRIBIN field:
   % Set BIN boundaries and make initial allocation of TRI_LIST (here for SEGS) for each bin:
   DELBIN = LXI / BODINT_PLANE.TBAXIS(AXIS).N_BINS;
   BODINT_PLANE.TBAXIS(AXIS).DELBIN = DELBIN;
   for IBIN=1:BODINT_PLANE.TBAXIS(AXIS).N_BINS
      BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).X_LOW  = BODINT_PLANE.BOX( LOW_IND,AXIS) + (IBIN-1)*DELBIN;
      BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).X_HIGH = BODINT_PLANE.BOX( LOW_IND,AXIS) + (IBIN  )*DELBIN;
      BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).NTL     = 0;
      BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).TRI_LIST= zeros(DELTA_SEGBIN);
   end
   % Finally, populate TRI_LIST (here for SEGS) for AXIS bins:
   for ISEG=1:BODINT_PLANE.NSEGS
      XIV(NOD1:NOD2) = BODINT_PLANE.XYZ(AXIS,BODINT_PLANE.SEGS(NOD1:NOD2,ISEG));
      XIV_LO  = min(XIV(NOD1:NOD2)); XIV_HI = max(XIV(NOD1:NOD2));
      AVAL   = (XIV_LO-GEOMEPS-BODINT_PLANE.BOX(LOW_IND,AXIS))/DELBIN;
      ILO_BIN= max(1,ceil(sign(AVAL)*min(2*BODINT_PLANE.TBAXIS(AXIS).N_BINS,abs(AVAL))));
      AVAL   = (XIV_HI+GEOMEPS-BODINT_PLANE.BOX(LOW_IND,AXIS))/DELBIN;
      IHI_BIN= min(BODINT_PLANE.TBAXIS(AXIS).N_BINS, ...
               ceil(sign(AVAL)*min(2*BODINT_PLANE.TBAXIS(AXIS).N_BINS,abs(AVAL))));
      for IBIN=ILO_BIN:IHI_BIN
         NTL = BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).NTL + 1;
         % Add Triangle index to BINs TRI_LIST
         BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).NTL = NTL;
         BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).TRI_LIST(NTL) = ISEG;
      end
   end
end

% Add Segments intersections:
for IBIN=1:BODINT_PLANE.TBAXIS(AXIS).N_BINS
   NTL = BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).NTL;
   % Now double loop, cost O(1/2*NTL^2):
   for BISEG=1:NTL
      ISEGV(EDG1) = BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).TRI_LIST(BISEG);
      
      SEGV(NOD1:NOD2,EDG1) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEGV(EDG1));
      P1(IAXIS:JAXIS)      = [ BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1,EDG1)), ...
                               BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1,EDG1))];
      D1(IAXIS:JAXIS)      = [ BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD2,EDG1)), ...
                               BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD2,EDG1))];
      D1 = D1 - P1;
      
      S1_X2_MIN=min(BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1:NOD2,EDG1)));
      S1_X2_MAX=max(BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1:NOD2,EDG1)));
      S1_X3_MIN=min(BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1:NOD2,EDG1)));
      S1_X3_MAX=max(BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1:NOD2,EDG1)));

      
      for BIISEG=BISEG+1:NTL
         ISEGV(EDG2) = BODINT_PLANE.TBAXIS(AXIS).TRIBIN(IBIN).TRI_LIST(BIISEG);
         SEGV(NOD1:NOD2,EDG2) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEGV(EDG2));
         P2(IAXIS:JAXIS)      = [ BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1,EDG2)), ...
                                  BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1,EDG2))];
         D2(IAXIS:JAXIS)      = [ BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD2,EDG2)), ...
                                  BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD2,EDG2))];
         D2 = D2 - P2;

         % Tests for quick discard:
         if( max(BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1:NOD2,EDG2)))+GEOMEPS < S1_X2_MIN); continue; end
         if( min(BODINT_PLANE.XYZ(X2AXIS,SEGV(NOD1:NOD2,EDG2)))-GEOMEPS > S1_X2_MAX); continue; end
         if( max(BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1:NOD2,EDG2)))+GEOMEPS < S1_X3_MIN); continue; end
         if( min(BODINT_PLANE.XYZ(X3AXIS,SEGV(NOD1:NOD2,EDG2)))-GEOMEPS > S1_X3_MAX); continue; end
         
         % Test for segment-segment intersection:
         [SVARV,SLENV,INT_FLG]=GET_SEGSEG_INTERSECTION(P1,D1,P2,D2);

         % Now discard repeated intersections:
         % If crossing is already defined in SEG don't add:
         for ICROSS=1:INT_FLG
           for ISX = EDG1:EDG2 
               SBOD = SVARV(ICROSS,ISX);           
               % Discard intersections already present in segment, including ends:
               INLIST=false;
               for ISVAR=1:BODINT_PLANE.NBCROSS(ISEGV(ISX))
                  if ( abs(SBOD-BODINT_PLANE.SVAR(ISVAR,ISEGV(ISX))) < GEOMEPS )
                      INLIST=true;
                      break
                  end
               end
               if (INLIST); continue; end
               
               % Add crossing to BODINT_PLANE, insertion sort:
               NBCROSS  = BODINT_PLANE.NBCROSS(ISEGV(ISX)) + 1;
               BODINT_PLANE.SVAR(NBCROSS,ISEGV(ISX)) = 1./GEOMEPS;
               for IBCR=1:NBCROSS
                  if ( SBOD < BODINT_PLANE.SVAR(IBCR,ISEGV(ISX)) ); break; end
               end

               % Here copy from the back (updated nbcross) to the ibcr location:
               for IDUM = NBCROSS:-1:IBCR+1
                  BODINT_PLANE.SVAR(IDUM,ISEGV(ISX)) = BODINT_PLANE.SVAR(IDUM-1,ISEGV(ISX));
               end
               BODINT_PLANE.SVAR(IBCR,ISEGV(ISX)) = SBOD;
               BODINT_PLANE.NBCROSS(ISEGV(ISX))   = NBCROSS;

               % Here we have an intersection inside a segment, note it in
               % FACERT:
               if (ISX==EDG1)
                   % X2AXIS, X3AXIS location of intersection:
                   XY(IAXIS:JAXIS) = P1(IAXIS:JAXIS) + SBOD*D1(IAXIS:JAXIS)/norm(D1(IAXIS:JAXIS));
               else
                   % X2AXIS, X3AXIS location of intersection:
                   XY(IAXIS:JAXIS) = P2(IAXIS:JAXIS) + SBOD*D2(IAXIS:JAXIS)/norm(D2(IAXIS:JAXIS));
               end
               XPOS = XY(IAXIS);
               if (X2NOC==0)
                   JJ2_LO = floor((XPOS-GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
                   JJ2_HI = floor((XPOS+GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
                   if (all([JJ2_LO JJ2_HI] < X2LO_CELL) || all([JJ2_LO JJ2_HI] > X2HI_CELL)); continue; end
                   JJ2_LO = max(JJ2_LO,X2LO_CELL); JJ2_HI = min(JJ2_HI,X2HI_CELL);
               else
                   FOUND_SEG = false;
                   JJ2_LO = -100;
                   JJ2_HI = -100;
                   for JJ2=X2LO_CELL:X2HI_CELL
                       % Check if XPOS is within this segment JJ2:
                       if ((XPOS-X2FACE(JJ2-1)) > -GEOMEPS && (X2FACE(JJ2)-XPOS) > -GEOMEPS)
                           if(JJ2_LO > -100)
                               JJ2_HI = JJ2;
                               break
                           else
                               JJ2_LO = JJ2;
                               JJ2_HI = JJ2;
                           end
                           FOUND_SEG=true;
                       end
                   end
                   if (~FOUND_SEG); continue; end
               end
               XPOS = XY(JAXIS);
               if (X3NOC==0)
                   KK2_LO  = floor((XPOS-GEOMEPS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
                   KK2_HI  = floor((XPOS+GEOMEPS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
                   if (all([KK2_LO KK2_HI] < X3LO_CELL) || all([KK2_LO KK2_HI]  > X3HI_CELL)); continue; end
                   KK2_LO = max(KK2_LO,X3LO_CELL); KK2_HI = min(KK2_HI,X3HI_CELL); 
               else
                   FOUND_SEG=false;
                   KK2_LO = -100;
                   KK2_HI = -100;
                   for KK2=X3LO_CELL:X3HI_CELL
                       % Check if XPOS is within this segment KK2:
                       if ((XPOS-X3FACE(KK2-1)) > -GEOMEPS && (X3FACE(KK2)-XPOS) > -GEOMEPS)
                           if(KK2_LO > -100)
                               KK2_HI = KK2;
                               break
                           else
                               KK2_LO = KK2;
                               KK2_HI = KK2;
                           end
                           FOUND_SEG=true;
                       end
                   end
                   if (~FOUND_SEG); continue; end
               end
               
               for KK2=KK2_LO:KK2_HI
                   for JJ2=JJ2_LO:JJ2_HI
                       FACERT(JJ2,KK2) = 1;
                   end
               end

           end
         end
         
      end
   end
end

% Loop nodes and test in SEG_NODES: if more than 2 segments end in the
% node, note it in FACERT.
MAX_SEG_NODE = max(SEGS_NODE(1:BODINT_PLANE.NNODS));
ISEG_NODE    = zeros(MAX_SEG_NODE+1,BODINT_PLANE.NNODS);
ANGS_NODE    = zeros(MAX_SEG_NODE  ,BODINT_PLANE.NNODS);
for ISEG=1:BODINT_PLANE.NSEGS
   % End nodes to cross:
   if( any(BODINT_PLANE.SEGTYPE(NOD1:NOD2,ISEG)~=IBM_GASPHASE) )
      SEG(NOD1:NOD2) = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
      DX2 = BODINT_PLANE.XYZ(X2AXIS,SEG(NOD2))-BODINT_PLANE.XYZ(X2AXIS,SEG(NOD1));
      DX3 = BODINT_PLANE.XYZ(X3AXIS,SEG(NOD2))-BODINT_PLANE.XYZ(X3AXIS,SEG(NOD1));
      for INOD=NOD1:NOD2
         % Compute angle, for NOD2 the seg andgle is -ANG.
         ANG=(NOD2-INOD)*atan2(DX3,DX2)+(INOD-NOD1)*atan2(-DX3,-DX2);
         if(ANG < 0.); ANG=ANG+2*pi; end % Make angle from 0 to 2*pi.        
         % Insert-add segment into ISEG_NODE depending on angle value:
         NSN                          = ISEG_NODE(1,SEG(INOD));
         ISEG_NODE(1      ,SEG(INOD)) = NSN+1;
         FOUND=false;
         ISEG2=1;
         if(NSN>0)
             for ISEG2=1:NSN
                 if(ANGS_NODE(ISEG2,SEG(INOD)) > ANG)
                     FOUND=true;
                     break
                 end
             end
         else
            ISEG_NODE(ISEG2+1,SEG(INOD)) =  ISEG;
            ANGS_NODE(ISEG2  ,SEG(INOD)) =   ANG;    
            continue
         end
         if(FOUND)
            for ISEG3=NSN+1:-1:ISEG2+1
                ISEG_NODE(ISEG3+1,SEG(INOD)) = ISEG_NODE(ISEG3  ,SEG(INOD));
                ANGS_NODE(ISEG3  ,SEG(INOD)) = ANGS_NODE(ISEG3-1,SEG(INOD));
            end
            ISEG_NODE(ISEG2+1,SEG(INOD)) =  ISEG;
            ANGS_NODE(ISEG2  ,SEG(INOD)) =   ANG;    
         else
            ISEG_NODE(ISEG2+2,SEG(INOD)) =  ISEG;
            ANGS_NODE(ISEG2+1,SEG(INOD)) =   ANG;    
         end
      end
   end
end

for INOD = 1:BODINT_PLANE.NNODS
    
   if(SEGS_NODE(INOD) < 3); continue; end
      
   % Test case of even number of segments:
   if (mod(SEGS_NODE(INOD),2)==0) % Case of even number of segments.
      % Test if circling around the node we have media discontinuity.
      NSN=ISEG_NODE(1,INOD);
      COUNT=0;
      for ISEG=ISEG_NODE(2:NSN+1,INOD)'
          COUNT=COUNT+1;
          SEG = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
          if (INOD==SEG(NOD2))
            CIRC_MED(COUNT) = BODINT_PLANE.SEGTYPE(NOD2,ISEG);
          else
            CIRC_MED(COUNT) = BODINT_PLANE.SEGTYPE(NOD1,ISEG);
          end
      end
      CIRC_MED(COUNT+1) =  CIRC_MED(1);
      CRS_FLG=false;
      for COUNT=1:NSN
          if(CIRC_MED(COUNT)==CIRC_MED(COUNT+1))
              CRS_FLG=true;
              break
          end
      end
      if(~CRS_FLG); continue; end
   end
             
   % X2AXIS, X3AXIS location of intersection:
   XY(IAXIS:JAXIS) = [BODINT_PLANE.XYZ(X2AXIS,INOD) BODINT_PLANE.XYZ(X3AXIS,INOD)]; 
    

   XPOS = XY(IAXIS);
   if (X2NOC==0)
       JJ2_LO = floor((XPOS-GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
       JJ2_HI = floor((XPOS+GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO)) + X2LO_CELL;
       if (all([JJ2_LO JJ2_HI] < X2LO_CELL) || all([JJ2_LO JJ2_HI] > X2HI_CELL)); continue; end
       JJ2_LO = max(JJ2_LO,X2LO_CELL); JJ2_HI = min(JJ2_HI,X2HI_CELL);
   else
       FOUND_SEG=false;
       JJ2_LO = -100;
       JJ2_HI = -100;
       for JJ2=X2LO_CELL:X2HI_CELL
           % Check if XPOS is within this segment JJ2:
           if ((XPOS-X2FACE(JJ2-1)) > -GEOMEPS && (X2FACE(JJ2)-XPOS) > -GEOMEPS)
               if(JJ2_LO > -100)
                   JJ2_HI = JJ2;
                   break
               else
                   JJ2_LO = JJ2;
                   JJ2_HI = JJ2;
               end
               FOUND_SEG=true;
           end
       end
       if (~FOUND_SEG); continue; end
   end
   XPOS = XY(JAXIS);
   if (X3NOC==0)
       KK2_LO  = floor((XPOS-GEOMEPS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
       KK2_HI  = floor((XPOS+GEOMEPS-X3FACE(X3LO))/DX3FACE(X3LO)) + X3LO_CELL;
       if (all([KK2_LO KK2_HI] < X3LO_CELL) || all([KK2_LO KK2_HI]  > X3HI_CELL)); continue; end
       KK2_LO = max(KK2_LO,X3LO_CELL); KK2_HI = min(KK2_HI,X3HI_CELL); 
   else
       FOUND_SEG=false;
       KK2_LO = -100;
       KK2_HI = -100;
       for KK2=X3LO_CELL:X3HI_CELL
           % Check if XPOS is within this segment KK2:
           if ((XPOS-X3FACE(KK2-1)) > -GEOMEPS && (X3FACE(KK2)-XPOS) > -GEOMEPS)
               if(KK2_LO > -100)
                   KK2_HI = KK2;
                   break
               else
                   KK2_LO = KK2;
                   KK2_HI = KK2;
               end
               FOUND_SEG=true;
           end
       end
       if (~FOUND_SEG); continue; end
   end

   for KK2=KK2_LO:KK2_HI
       for JJ2=JJ2_LO:JJ2_HI
           FACERT(JJ2,KK2) = 1;
       end
   end

end

% ISX=0;
% for KK2=X3LO_CELL:X3HI_CELL
%     for JJ2=X2LO_CELL:X2HI_CELL
%        if (FACERT(JJ2,KK2)) 
%            ISX = ISX+1;
%        end
%     end
% end
% disp(['X1AXIS,X1PLN,N FACERT=' num2str([X1AXIS,X1PLN,ISX])])

ierr=0;


return