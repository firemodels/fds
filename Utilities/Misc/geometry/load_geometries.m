function [GEOM,N_GEOM]=load_geometries(basedir,casename);

global MAX_DIM IAXIS JAXIS KAXIS NOD1 NOD2 NOD3 NOD4 LOW_IND HIGH_IND NODS_WSEL
global EDG1 EDG2 EDG3 GEOMEPS

FGEOM=load([basedir casename '_num_geometries.dat']);
N_GEOM=FGEOM(1);
GEOMEPS=FGEOM(2);
for IG=1:N_GEOM
   GEOM(IG).XYZ    =load([basedir casename '_geometry_' num2str(IG,'%4.4d') '_verts.dat']); 
   GEOM(IG).WSELEM =load([basedir casename '_geometry_' num2str(IG,'%4.4d') '_faces.dat']); 
   
   % Surf IDs:
   GEOM(IG).SURFS  = GEOM(IG).WSELEM(:,NOD4);
   GEOM(IG).WSELEM = GEOM(IG).WSELEM(:,NOD1:NOD3);
   
   % Verts:
   GEOM(IG).N_VERTS= length(GEOM(IG).XYZ(:,1));
   GEOM(IG).VERTS  = zeros(1,MAX_DIM*GEOM(IG).N_VERTS);
   COUNT=0;
   for IVERT=1:GEOM(IG).N_VERTS
       GEOM(IG).VERTS(COUNT+1:COUNT+3)=GEOM(IG).XYZ(IVERT,IAXIS:KAXIS);
       COUNT = COUNT+3;
   end
   
   % Faces:
   GEOM(IG).N_FACES= length(GEOM(IG).WSELEM(:,1));
   GEOM(IG).FACES  = zeros(1,NOD3*GEOM(IG).N_FACES);
   COUNT=0;
   for IFACE=1:GEOM(IG).N_FACES
       GEOM(IG).FACES(COUNT+1:COUNT+3)=GEOM(IG).WSELEM(IFACE,NOD1:NOD3);
       COUNT = COUNT+3;
   end   
   
   % FACECUBE
   GEOM(IG).FACECUBE = zeros(2,3,GEOM(IG).N_FACES);
   for IFACE=1:GEOM(IG).N_FACES
       XYZ=GEOM(IG).XYZ(GEOM(IG).WSELEM(IFACE,NOD1:NOD3),IAXIS:KAXIS);
       for AXIS=IAXIS:KAXIS
          GEOM(IG).FACECUBE(LOW_IND,AXIS,IFACE)  = min(XYZ(:,AXIS));
          GEOM(IG).FACECUBE(HIGH_IND,AXIS,IFACE) = max(XYZ(:,AXIS));
       end
   end
   
   % FACE_EDGES y EDGE_FACES:
   GEOM(IG).EDGES      = zeros(2,3*GEOM(IG).N_FACES);
   GEOM(IG).FACE_EDGES = zeros(3,GEOM(IG).N_FACES);
   GEOM(IG).EDGE_FACES = zeros(5,3*GEOM(IG).N_FACES);
   NWSEDG = 0;
   for IWSEL=1:GEOM(IG).N_FACES

      WSELEM(NOD1:NOD3) = GEOM(IG).FACES(NODS_WSEL*(IWSEL-1)+1:NODS_WSEL*IWSEL);

      if(mod(IWSEL,10000)==0) 
          disp(['1. IWSEL,N_FACES=' num2str([IWSEL GEOM(IG).N_FACES])])
      end
      
      for IEDGE=EDG1:EDG3 % For each face halfedge

         SEG(NOD1:NOD2) = WSELEM(NOD1:NOD2); % Get halfedge

         % Test triangles edge iedge is already on list
         % GEOM(IG).EDGES. Makes use of fact that two triangles
         % sharing an edge have opposite connectivity for it (right hand
         % rule for connectivity for normal outside solid).

         INLIST = false;
%          for IEDLIST=1:NWSEDG
%             % Check if opposite halfedge already in list.
%             if ( (SEG(NOD1) == GEOM(IG).EDGES(NOD2,IEDLIST)) && ...
%                  (SEG(NOD2) == GEOM(IG).EDGES(NOD1,IEDLIST)) )
%                INLIST = true;
%                break
%             end
%          end
         I1=find(GEOM(IG).EDGES(NOD2,1:NWSEDG)==SEG(NOD1));
         I2=find(GEOM(IG).EDGES(NOD1,1:NWSEDG)==SEG(NOD2));
         if (~isempty(I1) && ~isempty(I2))
             for II=1:length(I1)
                 II2=find(I2(:)==I1(II));
                 if(~isempty(II2))
                    INLIST = true;
                    IEDLIST= I1(II);
                    break
                 end
             end
         end
         
         if (INLIST) % Opposite halfedge already in list
            % Couple halfedge with its pair
            GEOM(IG).EDGE_FACES(1,IEDLIST)   = 2;
            GEOM(IG).EDGE_FACES(4,IEDLIST)   = IWSEL;
            GEOM(IG).EDGE_FACES(5,IEDLIST)   = IEDGE;
            GEOM(IG).FACE_EDGES(IEDGE,IWSEL) = IEDLIST;
         else % Opposite halfedge not in list, add a new entry
            NWSEDG = NWSEDG + 1;
            GEOM(IG).EDGES(NOD1:NOD2,NWSEDG) = SEG(NOD1:NOD2);
            GEOM(IG).EDGE_FACES(1,NWSEDG)    = 1;
            GEOM(IG).EDGE_FACES(2,NWSEDG)    = IWSEL;
            GEOM(IG).EDGE_FACES(3,NWSEDG)    = IEDGE;
            GEOM(IG).FACE_EDGES(IEDGE,IWSEL) = NWSEDG;
         end

         
         WSELEM=wshift('1D',WSELEM,1);
         
      end
   end

   
   
   
   GEOM(IG).N_EDGES = NWSEDG;
   GEOM(IG).FACES_NORMAL = zeros(MAX_DIM,GEOM(IG).N_FACES);
   GEOM(IG).FACES_AREA   = zeros(1,GEOM(IG).N_FACES);
   GEOM(IG).GEOM_VOLUME    = 0.;
   GEOM(IG).GEOM_AREA      = 0.;
   GEOM(IG).GEOM_XYZCEN(:) = 0.;

   % Compute normal, area and volume:
   SQAREA(IAXIS:KAXIS) = 0.;
   for IWSEL=1:GEOM(IG).N_FACES

      WSELEM(NOD1:NOD3) = GEOM(IG).FACES(NODS_WSEL*(IWSEL-1)+1:NODS_WSEL*IWSEL);

      if(mod(IWSEL,10000)==0) 
          disp(['2. IWSEL,N_FACES=' num2str([IWSEL GEOM(IG).N_FACES])])
      end      
      
      % Triangles NODES coordinates:
      for INOD=NOD1:NOD3
         XYZV(IAXIS:KAXIS,INOD) = GEOM(IG).VERTS(MAX_DIM*(WSELEM(INOD)-1)+1:MAX_DIM*WSELEM(INOD));
      end

      V12(IAXIS:KAXIS) = XYZV(IAXIS:KAXIS,NOD2) - XYZV(IAXIS:KAXIS,NOD1);
      V23(IAXIS:KAXIS) = XYZV(IAXIS:KAXIS,NOD3) - XYZV(IAXIS:KAXIS,NOD2);
      V31(IAXIS:KAXIS) = XYZV(IAXIS:KAXIS,NOD1) - XYZV(IAXIS:KAXIS,NOD3);

      % Cross V12 x V23:
      WSNORM(IAXIS) = V12(JAXIS)*V23(KAXIS) - V12(KAXIS)*V23(JAXIS);
      WSNORM(JAXIS) = V12(KAXIS)*V23(IAXIS) - V12(IAXIS)*V23(KAXIS);
      WSNORM(KAXIS) = V12(IAXIS)*V23(JAXIS) - V12(JAXIS)*V23(IAXIS);

      MGNRM = sqrt( WSNORM(IAXIS)^2. + WSNORM(JAXIS)^2. + WSNORM(KAXIS)^2. );

      XCEN  = (XYZV(IAXIS,NOD1) + XYZV(IAXIS,NOD2) + XYZV(IAXIS,NOD3)) / 3.;

      % Assign to GEOM:
      GEOM(IG).FACES_NORMAL(IAXIS:KAXIS,IWSEL) = WSNORM(IAXIS:KAXIS) * MGNRM^(-1.);
      GEOM(IG).FACES_AREA(IWSEL) = MGNRM/2.;
     
      % Total Area and Volume for GEOMETRY(IG).
      GEOM(IG).GEOM_AREA  = GEOM(IG).GEOM_AREA  + GEOM(IG).FACES_AREA(IWSEL);
      GEOM(IG).GEOM_VOLUME= GEOM(IG).GEOM_VOLUME+ ... % Divergence theorem with F = x i, assumes we have a volume.
      GEOM(IG).FACES_NORMAL(IAXIS,IWSEL)*XCEN*GEOM(IG).FACES_AREA(IWSEL);

      % Define Centroid:
      X12(IAXIS:KAXIS) = 0.5*(XYZV(IAXIS:KAXIS,NOD1) + XYZV(IAXIS:KAXIS,NOD2))';
      X23(IAXIS:KAXIS) = 0.5*(XYZV(IAXIS:KAXIS,NOD2) + XYZV(IAXIS:KAXIS,NOD3))';
      X31(IAXIS:KAXIS) = 0.5*(XYZV(IAXIS:KAXIS,NOD3) + XYZV(IAXIS:KAXIS,NOD1))';
      
      % dot(i,nc) int(x^2)dA, dot(j,nc) int(y^2)dA, dot(k,nc) int(z^2)dA
      for IX=IAXIS:KAXIS
         INT2 = (X12(IX)^2. + X23(IX)^2. + X31(IX)^2.) / 3.;
         SQAREA(IX) = SQAREA(IX) + GEOM(IG).FACES_NORMAL(IX,IWSEL)*INT2*GEOM(IG).FACES_AREA(IWSEL); % Midpt rule.
      end

   
   end 
   
   % Geometry Centroid:
   for IX=IAXIS:KAXIS
      GEOM(IG).GEOM_XYZCEN(IX) = SQAREA(IX) / (2. * GEOM(IG).GEOM_VOLUME);
   end
   
    % Check we don't have repeated nodes:
    for IVERT=1:GEOM(IG).N_VERTS
        DIFF=abs(GEOM(IG).XYZ(:,IAXIS)-GEOM(IG).XYZ(IVERT,IAXIS)) + ...
             abs(GEOM(IG).XYZ(:,JAXIS)-GEOM(IG).XYZ(IVERT,JAXIS)) + ...
             abs(GEOM(IG).XYZ(:,KAXIS)-GEOM(IG).XYZ(IVERT,KAXIS));
        [vdf,idf]=find(DIFF < GEOMEPS);
        if (length(vdf)>1)
            disp(['Found same nodes: ' num2str([IVERT idf])])
        end
    end

    % Check we don't have repeated faces:
    for IFACE=1:GEOM(IG).N_FACES
        DIFF=abs(GEOM(IG).WSELEM(:,NOD1)-GEOM(IG).WSELEM(IFACE,NOD1)) + ...
             abs(GEOM(IG).WSELEM(:,NOD2)-GEOM(IG).WSELEM(IFACE,NOD2)) + ...
             abs(GEOM(IG).WSELEM(:,NOD3)-GEOM(IG).WSELEM(IFACE,NOD3));
        [vdf,idf]=find(DIFF < GEOMEPS);
        if (length(vdf)>1)
            disp(['Found same faces: ' num2str([IFACE idf])])
        end
    end
    
    % Check that all edges attached to 2 faces:
    for IEDGE=1:GEOM(IG).N_EDGES 
        if(GEOM(IG).EDGE_FACES(1,IEDGE) < 2)
            disp(['IEDGE=' num2str(IEDGE) ' has < 2 faces.'])
            pause
        end
    end
         
end

return