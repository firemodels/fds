function [ierr]=GET_CARTEDGE_CUTEDGES(X1AXIS,X2AXIS,X3AXIS,XIAXIS,XJAXIS,XKAXIS, ...
                                      NM,X2LO_CELL,X2HI_CELL,INDX1,KK)
                                 
global MAX_DIM IAXIS JAXIS KAXIS FCELL NOD1 NOD2 GEOMEPS MESHES
global IBM_SOLID IBM_CUTCFE IBM_GASPHASE IBM_UNDEFINED IBM_VGSC IBM_EGSC IBM_ECRS IBM_IDCE
global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS
global X2NOC IBM_SS IBM_GS IBM_SG
global X2CELL DX2CELL X1FACE X2FACE X3FACE

ierr=1;
INDXI = zeros(1,KAXIS);
% Set initially edges with MESHES(NM).VERTVAR vertices == IBM_SOLID to IBM_SOLID status:
for JJ=X2LO_CELL:X2HI_CELL
    % Vert at index JJ-FCELL:
    INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ-FCELL, KK ]; % Local x1,x2,x3
    INDI=INDXI(XIAXIS);
    INDJ=INDXI(XJAXIS);
    INDK=INDXI(XKAXIS);
    
    % Vert at index JJ-FCELL+1:
    INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ-FCELL+1, KK ]; % Local x1,x2,x3
    INDI1=INDXI(XIAXIS);
    INDJ1=INDXI(XJAXIS);
    INDK1=INDXI(XKAXIS);

    % Edge at index JJ:
    INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ, KK ]; % Local x1,x2,x3
    INDIE=INDXI(XIAXIS);
    INDJE=INDXI(XJAXIS);
    INDKE=INDXI(XKAXIS);

    if ((MESHES(NM).VERTVAR(INDI ,INDJ ,INDK ,IBM_VGSC) == IBM_SOLID) && ...
        (MESHES(NM).VERTVAR(INDI1,INDJ1,INDK1,IBM_VGSC) == IBM_SOLID) )     
         MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_EGSC,X2AXIS) = IBM_SOLID;
    end
end

NEDGECROSS_OLD = MESHES(NM).N_EDGE_CROSS;

% Edges with Crossings not on VERTICES:
for ICRS=1:IBM_N_CRS

    % Skip SOLID-SOLID intersections, as there is no media crossing:
    if (IBM_IS_CRS(ICRS) == IBM_SS); continue; end

    % Check location on grid of crossing:
    % See if crossing is exactly on a Cartesian cell vertex:
    if (X2NOC==0)
       % Optimized for UG:
       JSTR = floor( (IBM_SVAR_CRS(ICRS)-GEOMEPS-X2CELL(X2LO_CELL))/DX2CELL(X2LO_CELL) ) + X2LO_CELL;
       % Discard cut-edges on Cartesian edges laying > X2HI_CELL.
       if (JSTR < X2LO_CELL-1); continue; end
       if (JSTR > X2HI_CELL+1); continue; end

       JJ    = JSTR;
       DELJJ = abs(X2CELL(JJ)-IBM_SVAR_CRS(ICRS)) - DX2CELL(X2LO_CELL)/2.;
       % Crossing on Vertex?
       if ( abs(DELJJ) < GEOMEPS ) % Add crossing to two edges:
           JJLOW=0; JJHIGH=1;
       elseif ( DELJJ < -GEOMEPS ) % Crossing in jj Edge.
           JJLOW=0; JJHIGH=0;
       elseif ( DELJJ > GEOMEPS )  % Crossing in jj+1 Edge.
           JJLOW=1; JJHIGH=1;
       end
    else
       FOUND_EDGE=false;
       JJLOW = -1000000;
       JJHIGH=  1000000;
       for JJ=X2LO_CELL-1:X2HI_CELL
          DELJJ = IBM_SVAR_CRS(ICRS)-X2CELL(JJ);
          XVJJ  = X2CELL(JJ) + DX2CELL(JJ)/2.;
          DELJJ1= IBM_SVAR_CRS(ICRS)-X2CELL(JJ+1);
          % First two edges:
          if (abs(IBM_SVAR_CRS(ICRS)-XVJJ) < GEOMEPS) % Both JJ and JJ+1
             FOUND_EDGE=true;
             JJLOW=0; JJHIGH=1;
             break
          elseif (abs(DELJJ) <   DX2CELL(JJ)/2.) % JJ
             FOUND_EDGE=true;
             JJLOW=0; JJHIGH=0;
             break
          elseif (abs(DELJJ1)< DX2CELL(JJ+1)/2.) % JJ+1
             FOUND_EDGE=true;
             JJLOW=1; JJHIGH=1;
             break
          end
       end
       if (~FOUND_EDGE); continue; end
    end

    for JJADD=JJLOW:JJHIGH
        % Edge in the left:
        % Edge at index JJ or JJ+1:
        INDXI(IAXIS:KAXIS) = [ INDX1(X1AXIS), JJ+JJADD, KK ]; % Local x1,x2,x3
        INDIE=INDXI(XIAXIS);
        INDJE=INDXI(XJAXIS);
        INDKE=INDXI(XKAXIS);

        % Set MESHES(NM).ECVAR(IE,JE,KE,IBM_EGSC,X2AXIS) = IBM_CUTCFE:
        ICROSS = MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_ECRS,X2AXIS);

        if ( ICROSS > 0 ) % Edge has crossings already.

            % Populate EDGECROSS struct:
            NCROSS = MESHES(NM).EDGE_CROSS(ICROSS).NCROSS + 1;
            MESHES(NM).EDGE_CROSS(ICROSS).NCROSS       = NCROSS;
            MESHES(NM).EDGE_CROSS(ICROSS).SVAR(NCROSS) = IBM_SVAR_CRS(ICRS);
            MESHES(NM).EDGE_CROSS(ICROSS).ISVAR(NCROSS)= IBM_IS_CRS(ICRS);

        else % No crossings yet.

            NEDGECROSS = MESHES(NM).N_EDGE_CROSS + 1;
            MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_EGSC,X2AXIS) = IBM_CUTCFE;
            MESHES(NM).N_EDGE_CROSS                      = NEDGECROSS;
            MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_ECRS,X2AXIS) = NEDGECROSS;

            % Populate EDGECROSS struct:
            NCROSS = 1;
            MESHES(NM).EDGE_CROSS(NEDGECROSS).NCROSS       = NCROSS;
            MESHES(NM).EDGE_CROSS(NEDGECROSS).SVAR(NCROSS) = IBM_SVAR_CRS(ICRS);
            MESHES(NM).EDGE_CROSS(NEDGECROSS).ISVAR(NCROSS)= IBM_IS_CRS(ICRS);
            MESHES(NM).EDGE_CROSS(NEDGECROSS).IJK(1:4) = [ INDIE, INDJE, INDKE, X2AXIS ];

        end

    end

end

% Now Define MESHES(NM).CUT_EDGE for IBM_GASPHASE cut-edges:
% First: Run over crossings and set MESHES(NM).IBM_CUT_EDGES:
for ICROSS=NEDGECROSS_OLD+1:MESHES(NM).N_EDGE_CROSS

   % Discard edge outside of blocks ranges for ray on x2axis:
   if ( (MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS) < X2LO_CELL) || ...
        (MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS) > X2HI_CELL) ); continue; end

   NCROSS = MESHES(NM).EDGE_CROSS(ICROSS).NCROSS;

   % Edge Location in x1,x2,x3 axes:
   % Vert at index JJ-FCELL:
   INDXI(IAXIS:KAXIS) = [ MESHES(NM).EDGE_CROSS(ICROSS).IJK(X1AXIS),       ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL, ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X3AXIS) ]; % Local x1,x2,x3
   INDI=INDXI(XIAXIS);
   INDJ=INDXI(XJAXIS);
   INDK=INDXI(XKAXIS);
   % Vert at index JJ-FCELL+1:
   INDXI(IAXIS:KAXIS) = [ MESHES(NM).EDGE_CROSS(ICROSS).IJK(X1AXIS),           ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL+1,   ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X3AXIS) ]; % Local x1,x2,x3
   INDI1=INDXI(XIAXIS);
   INDJ1=INDXI(XJAXIS);
   INDK1=INDXI(XKAXIS);
   % Edge at index jj:
   INDXI(IAXIS:KAXIS) = [ MESHES(NM).EDGE_CROSS(ICROSS).IJK(X1AXIS),    ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS),    ...
                          MESHES(NM).EDGE_CROSS(ICROSS).IJK(X3AXIS) ]; % Local x1,x2,x3
   INDIE=INDXI(XIAXIS); % i.e. MESHES(NM).EDGE_CROSS(ICROSS).IJK(IAXIS), etc.
   INDJE=INDXI(XJAXIS);
   INDKE=INDXI(XKAXIS);

   % Discard Edge with one EDGECROSS and both vertices having VERTVAR = IBM_SOLID:
   % The crossing is on one of the edge vertices.
   if ( (NCROSS == 1)                                                 && ...
        (MESHES(NM).VERTVAR(INDI ,INDJ ,INDK ,IBM_VGSC) == IBM_SOLID) && ...
        (MESHES(NM).VERTVAR(INDI1,INDJ1,INDK1,IBM_VGSC) == IBM_SOLID) )

      MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_EGSC,X2AXIS) = IBM_SOLID;
      continue

   end

   % Discard cases for edge with two crossings:
   if ( NCROSS == 2 )

      VSOLID = (MESHES(NM).VERTVAR(INDI ,INDJ ,INDK ,IBM_VGSC) == IBM_SOLID) && ...
               (MESHES(NM).VERTVAR(INDI1,INDJ1,INDK1,IBM_VGSC) == IBM_SOLID);

      % Test if crossings lay on same location + solid vertices:
      DIF  = ( MESHES(NM).EDGE_CROSS(ICROSS).SVAR(2) - ...
               MESHES(NM).EDGE_CROSS(ICROSS).SVAR(1) ) < GEOMEPS;
      if (DIF && VSOLID)
         MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_EGSC,X2AXIS) = IBM_SOLID;
         continue
      end

      DIF  = (abs(X2FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL  ) -     ...
                         MESHES(NM).EDGE_CROSS(ICROSS).SVAR(1)) < GEOMEPS)  &&    ...
             (abs(X2FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL+1) -     ...
                         MESHES(NM).EDGE_CROSS(ICROSS).SVAR(2)) < GEOMEPS);

      VFLUID  = (MESHES(NM).EDGE_CROSS(ICROSS).ISVAR(1) == IBM_GS)  && ...
                (MESHES(NM).EDGE_CROSS(ICROSS).ISVAR(2) == IBM_SG);

      if (DIF && VSOLID && VFLUID)
         MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_EGSC,X2AXIS) = IBM_SOLID;
         continue
      end

   end

   % New CUT_EDGE struct for this edge:
   NCUTEDGE = MESHES(NM).N_CUTEDGE_MESH + 1;
   MESHES(NM).N_CUTEDGE_MESH                          = NCUTEDGE;
   MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_IDCE,X2AXIS)= NCUTEDGE;

   MESHES(NM).CUT_EDGE(NCUTEDGE).STATUS           = IBM_GASPHASE;
   MESHES(NM).CUT_EDGE(NCUTEDGE).IJK(1:MAX_DIM+1) = MESHES(NM).EDGE_CROSS(ICROSS).IJK(1:MAX_DIM+1);
   MESHES(NM).CUT_EDGE(NCUTEDGE).IJK(MAX_DIM+2)   = IBM_UNDEFINED;    % No need to define type of CUT_EDGE
                                                                      % (is IBM_GASPHASE).
   % First Vertices:
   NVERT = NCROSS + 2;
   MESHES(NM).CUT_EDGE(NCUTEDGE).NVERT = NVERT;
   X123VERT(IAXIS:KAXIS,1:NVERT) = 0.;
   X123VERT(IAXIS,1:NVERT)   = X1FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X1AXIS));
   X123VERT(JAXIS,1)         = X2FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL);
   X123VERT(JAXIS,2:NCROSS+1)= MESHES(NM).EDGE_CROSS(ICROSS).SVAR(1:NCROSS);
   X123VERT(JAXIS,NVERT)     = X2FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X2AXIS)-FCELL+1);
   X123VERT(KAXIS,1:NVERT)   = X3FACE(MESHES(NM).EDGE_CROSS(ICROSS).IJK(X3AXIS));

   % Allocate new edge XYZVERT, CEELEM, INDSEG:
   for IVERT=1:MESHES(NM).CUT_EDGE(NCUTEDGE).NVERT
      MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(IAXIS:KAXIS,IVERT) = ...
                        X123VERT( [ XIAXIS, XJAXIS, XKAXIS ] ,IVERT);
   end

   % Now Cut Edges:
   % This assumes crossings are ordered for increasing svar, no repeated svar:
   NEDGE = 0;
   MESHES(NM).CUT_EDGE(NCUTEDGE).NEDGE = NEDGE;
   for IVERT=1:MESHES(NM).CUT_EDGE(NCUTEDGE).NVERT-1

      % Drop zero length edge (in x2 local dir):
      if (abs(X123VERT(JAXIS,IVERT)-X123VERT(JAXIS,IVERT+1)) < GEOMEPS); continue; end

      % Define if the cut-edge is gasphase:
      % Ray tracing for the center of the cut-edge most robust.
      XCEN  = 0.5*(MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(IAXIS,IVERT  ) + ...
                   MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(IAXIS,IVERT+1));
      YCEN  = 0.5*(MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(JAXIS,IVERT  ) + ...
                   MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(JAXIS,IVERT+1));
      ZCEN  = 0.5*(MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(KAXIS,IVERT  ) + ...
                   MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT(KAXIS,IVERT+1));
      XYZCEN(IAXIS:KAXIS) = [ XCEN, YCEN, ZCEN ];

      % Do a SOLID crossing count up to XYZcen(x2axis):
      SCEN=XYZCEN(X2AXIS);
      [IS_GASPHASE]=GET_IS_GASPHASE(SCEN);

      if ( IS_GASPHASE )
         NEDGE = NEDGE + 1;
         % Test for size of CEELEM, INDSEG, if smaller than NEDGE reallocate:
         MESHES(NM).CUT_EDGE(NCUTEDGE).NEDGE = NEDGE;
         MESHES(NM).CUT_EDGE(NCUTEDGE).CEELEM(NOD1:NOD2,NEDGE) = [ IVERT, IVERT+1 ];
      end
   end

   if (MESHES(NM).CUT_EDGE(NCUTEDGE).NEDGE == 0) % REWIND
      MESHES(NM).CUT_EDGE(NCUTEDGE).XYZVERT = 0.;
      MESHES(NM).CUT_EDGE(NCUTEDGE).CEELEM  = IBM_UNDEFINED;
      MESHES(NM).CUT_EDGE(NCUTEDGE).INDSEG  = IBM_UNDEFINED;
      NCUTEDGE = NCUTEDGE - 1;
      MESHES(NM).N_CUTEDGE_MESH             = NCUTEDGE;
      MESHES(NM).ECVAR(INDIE,INDJE,INDKE,IBM_IDCE,X2AXIS) = 0;
   end

end

ierr=0;

return