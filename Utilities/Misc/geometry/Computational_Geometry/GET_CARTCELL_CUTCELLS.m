function [ierr]=GET_CARTCELL_CUTCELLS(NM)

global LOW_IND HIGH_IND NGUARD CCGUARD FCELL
global IAXIS JAXIS KAXIS MAX_DIM
global ILO_CELL IHI_CELL JLO_CELL JHI_CELL KLO_CELL KHI_CELL 
global XFACE YFACE ZFACE
global IBM_CGSC IBM_FGSC IBM_CUTCFE IBM_UNDEFINED IBM_GASPHASE IBM_SOLID IBM_IDCF IBM_IDCC
global NOD1 NOD2 NOD3 NOD4 
global MESHES  GEOM
global IBM_MAXVERTS_CELL IBM_NPARAM_CCFACE
global DELTA_VERT DELTA_EDGE DELTA_FACE DELTA_CELL
global IBM_FTYPE_RGGAS IBM_FTYPE_CFGAS IBM_FTYPE_CFINB
global XCELL YCELL ZCELL DXCELL DYCELL DZCELL
global GEOMEPS

global CELLRT

ierr=1;

SIZE_CEELEM_EDGFAC = DELTA_EDGE;
SIZE_CFELEM_EDGFAC = DELTA_FACE;
SIZE_CEELEM_FACEDG = DELTA_EDGE;
SIZE_CFELEM_FACEDG = DELTA_FACE;
SIZE_VERTS_FC      = DELTA_VERT;
SIZE_CFELEM_FC     = DELTA_FACE;
SIZE_FACE_CCELEM   = DELTA_FACE;
SIZE_CELL_CCELEM   = DELTA_CELL;

% Definition of cut-cells:
% For each cartesian cell being cut into one or several cut-cells (NCELL), fill
% entries on a MESHES(NM).CUT_CELL struct. On each local entry ICC:
% - Add number of faces that are boundary of cut-cell.
%    MESHES(NM).CUT_CELL(ICELL)%CCELEM(1:NFACE_CELL+1,ICC), ICC=1,...,MESHES(NM).CUT_CELL(ICELL)%NCELL
% - Add list of corresponding regular faces, or cut-faces in CUT_FACE:
%    + 5 Indexes:
%    MESHES(NM).CUT_CELL(ICELL)%FACES_LIST = [ FACE_TYPE      LOW/HIGH      AXIS       cei      icf   ]
%                                                  where in MESHES(NM).CUT_FACE(CEI), which icf.
% - Compute Volume properties for each disjoint volume, add an unknown
%   number for scalars, pressure, etc.

for IBNDINT=LOW_IND:HIGH_IND  % 1 refers to blocks internal cells, 2 refers to block guard cells.

switch(IBNDINT)
   case(LOW_IND)
   IJK_COUNT = zeros(IHI_CELL+NGUARD,JHI_CELL+NGUARD,KHI_CELL+NGUARD);   
   ILO = ILO_CELL; IHI = IHI_CELL;
   JLO = JLO_CELL; JHI = JHI_CELL;
   KLO = KLO_CELL; KHI = KHI_CELL;
   case(HIGH_IND)
   ILO = ILO_CELL-CCGUARD; IHI = IHI_CELL+CCGUARD;
   JLO = JLO_CELL-CCGUARD; JHI = JHI_CELL+CCGUARD;
   KLO = KLO_CELL-CCGUARD; KHI = KHI_CELL+CCGUARD;
end

% Loop on Cartesian cells, define cut cells and solid cells IBM_CGSC:
for K=KLO:KHI
   for J=JLO:JHI
      for I=ILO:IHI

         if( MESHES(NM).CCVAR(I,J,K,IBM_CGSC) ~= IBM_CUTCFE ); continue; end

         if( IJK_COUNT(I,J,K) ); continue; end; IJK_COUNT(I,J,K) = true;

         % Start with Cartesian Faces:
         % Face type of bounding Cartesian faces:
         FSID_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_FGSC,IAXIS);
         FSID_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_FGSC,JAXIS);
         FSID_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_FGSC,JAXIS);
         FSID_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_FGSC,KAXIS);
         FSID_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_FGSC,KAXIS);

         % Cut-face number of bounding Cartesian faces:
         IDCF_XYZ(LOW_IND ,IAXIS) = MESHES(NM).FCVAR(I-FCELL  ,J,K,IBM_IDCF,IAXIS);
         IDCF_XYZ(HIGH_IND,IAXIS) = MESHES(NM).FCVAR(I-FCELL+1,J,K,IBM_IDCF,IAXIS);
         IDCF_XYZ(LOW_IND ,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL  ,K,IBM_IDCF,JAXIS);
         IDCF_XYZ(HIGH_IND,JAXIS) = MESHES(NM).FCVAR(I,J-FCELL+1,K,IBM_IDCF,JAXIS);
         IDCF_XYZ(LOW_IND ,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL  ,IBM_IDCF,KAXIS);
         IDCF_XYZ(HIGH_IND,KAXIS) = MESHES(NM).FCVAR(I,J,K-FCELL+1,IBM_IDCF,KAXIS);

         % Local variables:
         % Geometric entities related to the Cartesian cell:
         NVERT_CELL = 0;
         NSEG_CELL  = 0;
         NFACE_CELL = 0;
         SEG_CELL   = IBM_UNDEFINED*ones(NOD2,SIZE_CEELEM_EDGFAC);
         FACE_CELL  = IBM_UNDEFINED*ones(SIZE_VERTS_FC,SIZE_CFELEM_FC);
         FACE_LIST  = IBM_UNDEFINED*ones(IBM_NPARAM_CCFACE,SIZE_CFELEM_FC);
         XYZVERT    = zeros(KAXIS,IBM_MAXVERTS_CELL);
         AREAVARS   = zeros(MAX_DIM+1,SIZE_CFELEM_FC);

         % Add Cartesian Regular faces + GASPHASE cut-faces + vertices:
         IED = I-FCELL; JED = J-FCELL; KED = K-FCELL;
         for MYAXIS=IAXIS:KAXIS
            switch(MYAXIS)
            case(IAXIS)

               XYZLH(IAXIS:KAXIS,NOD1,LOW_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD2,LOW_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD3,LOW_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD4,LOW_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';

               XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND)= [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND)= [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND)= [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND)= [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';

               AREAI    = DYCELL(J) * DZCELL(K);
               AREAVARSI(1:MAX_DIM+1,LOW_IND) =[-XFACE(IED  )*AREAI, -XFACE(IED  )^2.*AREAI, 0., 0. ]';
               AREAVARSI(1:MAX_DIM+1,HIGH_IND)=[ XFACE(IED+1)*AREAI,  XFACE(IED+1)^2.*AREAI, 0., 0. ]';
            case(JAXIS)

               XYZLH(IAXIS:KAXIS,NOD1,LOW_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD2,LOW_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD3,LOW_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD4,LOW_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';

               XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND)= [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND)= [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND)= [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND)= [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';

               AREAI    = DXCELL(I) * DZCELL(K);
               AREAVARSI(1:MAX_DIM+1,LOW_IND) =[ 0., 0., -YFACE(JED  )^2.*AREAI, 0. ]';
               AREAVARSI(1:MAX_DIM+1,HIGH_IND)=[ 0., 0.,  YFACE(JED+1)^2.*AREAI, 0. ]';
            case(KAXIS)

               XYZLH(IAXIS:KAXIS,NOD1,LOW_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD2,LOW_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD3,LOW_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';
               XYZLH(IAXIS:KAXIS,NOD4,LOW_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';

               XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND)= [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND)= [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND)= [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
               XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND)= [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';

               AREAI    = DXCELL(I) * DYCELL(J);
               AREAVARSI(1:MAX_DIM+1,LOW_IND) =[ 0., 0., 0., -ZFACE(KED  )^2.*AREAI ]';
               AREAVARSI(1:MAX_DIM+1,HIGH_IND)=[ 0., 0., 0.,  ZFACE(KED+1)^2.*AREAI ]';
            end

            CEI_AXIS(LOW_IND:HIGH_IND) = IDCF_XYZ(LOW_IND:HIGH_IND,MYAXIS)';

            for SIDE=LOW_IND:HIGH_IND
               % Low High face:
               if ( FSID_XYZ(SIDE,MYAXIS) == IBM_GASPHASE )

                  % Regular Face, build 4 vertices + face:
                  NP = 0;
                  NFACE_CELL = NFACE_CELL + 1;

                  % Here, reallocate FACE_LIST, AREAVARS, FACE_CELL if NFACE_CELL > SIZE_CFELEM_FC:
                  % Also no need to reallocate FACE_CELL vert dimension, as for regular cells vert size = 5.
                  % CALL REALLOCATE_LOCAL_FC_VARS
                  FACE_LIST(1:IBM_NPARAM_CCFACE,NFACE_CELL) = [ IBM_FTYPE_RGGAS, SIDE, MYAXIS, 0, 0 ]';
                  % IBM_FTYPE_RGGAS=0, regular face.
                  AREAVARS(1:MAX_DIM+1,NFACE_CELL) = AREAVARSI(1:MAX_DIM+1,SIDE);

                  % Vertices arranged normal out of cartesian cell:
                  for IP=NOD1:NOD4
                     % xl,yl,zl
                     XYZ(IAXIS:KAXIS) = XYZLH(IAXIS:KAXIS,IP,SIDE)';
                     [NVERT_CELL,INOD,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_CELL,XYZ,NVERT_CELL,XYZVERT);

                     NP = NP + 1;
                     FACE_CELL(1,NFACE_CELL)    =   NP;
                     FACE_CELL(NP+1,NFACE_CELL) = INOD;
                  end

               elseif (FSID_XYZ(SIDE,MYAXIS) == IBM_CUTCFE )

                  FCT = 2*SIDE-3; % 2*(side-3/2);
                  % GasPhase CUT_FACE, add all cut-faces on these Cartesian cell + nodes:
                  CEI = CEI_AXIS(SIDE);
                  for ICF=1:MESHES(NM).CUT_FACE(CEI).NFACE
                     NFACE_CELL = NFACE_CELL + 1;
                     % Here, reallocate FACE_LIST, AREAVARS, FACE_CELL if NFACE_CELL > SIZE_CFELEM_FC:
                     % CALL REALLOCATE_LOCAL_FC_VARS
                     % Also reallocate FACE_CELL vert dimension, if needed.
                     NP = MESHES(NM).CUT_FACE(CEI).CFELEM(1,ICF);
                     % CALL REALLOCATE_FACE_CELL_VERTS

                     FACE_LIST(1:IBM_NPARAM_CCFACE,NFACE_CELL) = [ IBM_FTYPE_CFGAS, SIDE, MYAXIS, CEI, ICF ]';
                     % IBM_FTYPE_CFGAS=1
                     AREAVARS(1:MAX_DIM+1,NFACE_CELL) =[ MESHES(NM).CUT_FACE(CEI).INXAREA(ICF),   ...
                                                         MESHES(NM).CUT_FACE(CEI).INXSQAREA(ICF), ...
                                                         MESHES(NM).CUT_FACE(CEI).JNYSQAREA(ICF), ...
                                                         MESHES(NM).CUT_FACE(CEI).KNZSQAREA(ICF) ]'*FCT;
                                                         % FCT considers Normal out.
                     FACE_CELL(1,NFACE_CELL) = NP;
                     for IP=2:NP+1
                        FNOD             = MESHES(NM).CUT_FACE(CEI).CFELEM(IP,ICF);
                        XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_FACE(CEI).XYZVERT(IAXIS:KAXIS,FNOD)';
                        [NVERT_CELL,INOD,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_CELL,XYZ,NVERT_CELL,XYZVERT);
                        FACE_CELL(IP,NFACE_CELL) = INOD;
                     end
                  end
               end
            end
         end

         N_GAS_CFACES = NFACE_CELL;

         % Now add INBOUNDARY faces of the cell:
         CEI = MESHES(NM).CCVAR(I,J,K,IBM_IDCF);
         if ( CEI > 0 )
            FCT = -1.;
            for ICF=1:MESHES(NM).CUT_FACE(CEI).NFACE
               NFACE_CELL = NFACE_CELL + 1;
               % Here, reallocate FACE_LIST, AREAVARS, FACE_CELL if NFACE_CELL > SIZE_CFELEM_FC:
               % CALL REALLOCATE_LOCAL_FC_VARS
               % Also reallocate FACE_CELL, FACE_CELL_DUM vert dimension, if needed.
               NP = MESHES(NM).CUT_FACE(CEI).CFELEM(1,ICF);
               % CALL REALLOCATE_FACE_CELL_VERTS
               FACE_LIST(1:IBM_NPARAM_CCFACE,NFACE_CELL) = [ IBM_FTYPE_CFINB, 0, 0, CEI, ICF ]';
               % IBM_FTYPE_CFINB in Cart-cell.
               AREAVARS(1:MAX_DIM+1,NFACE_CELL) = [ MESHES(NM).CUT_FACE(CEI).INXAREA(ICF),   ...
                                                    MESHES(NM).CUT_FACE(CEI).INXSQAREA(ICF), ...
                                                    MESHES(NM).CUT_FACE(CEI).JNYSQAREA(ICF), ...
                                                    MESHES(NM).CUT_FACE(CEI).KNZSQAREA(ICF) ]'*FCT;
                                                    % Normal out of cut-cell.
               FACE_CELL(1,NFACE_CELL) = NP;
               for IP=2:NP+1
                  FNOD             = MESHES(NM).CUT_FACE(CEI).CFELEM(IP,ICF);
                  XYZ(IAXIS:KAXIS) = MESHES(NM).CUT_FACE(CEI).XYZVERT(IAXIS:KAXIS,FNOD)';
                  [NVERT_CELL,INOD,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_CELL,XYZ,NVERT_CELL,XYZVERT);
                  FACE_CELL(IP,NFACE_CELL) = INOD;
               end
               % At this point the face in face cell is ordered
               % throught the normal outside the body. Reorganize
               % to normal outside cut-cell (inside body).
               FACE_CELL_DUM(1:NP+1) = FACE_CELL(1:NP+1,NFACE_CELL)';
               for IP=2:NP+1
                  FACE_CELL(IP,NFACE_CELL) = FACE_CELL_DUM( (NP+1)+2-IP );
               end
            end
         end

         
%          if(I==36 && J==17 && K==24)
%             figure
%             hold on
%             for JCF=1:NFACE_CELL
%                NELEM  = FACE_CELL(1,JCF);
%                CFELEM = FACE_CELL(2:NELEM+1,JCF);
%            
%                if(FACE_LIST(1,JCF)==IBM_FTYPE_CFINB)
%                [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'r');
%                else
%                [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'b','Marker','o');
%                end
%                set(hp,'FaceAlpha',0.3) %,'EdgeAlpha',0.4) % .4    
%             end 
%             a=0.0005;
%             for IVERT=1:NVERT_CELL
%                 text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
%                      XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
%                 
%             end
%             axis equal; box on;
%             view([45 45])
%             xlabel('X')
%             ylabel('Y')
%             zlabel('Z')
%  
%             XYZVERT(:,[20 22])
%             pause
%          end
         
         

         % Here we have in XYZvert all the vertices that define the
         % cut-cells within Cartesian cell I,J,K. We have the faces,
         % boundary of said cut-cells in face_cell.
         % We have in face_list the list of cut-cell boundary faces
         % and if they are regular or cut-face.
         % We want to reorder face list, such that we have the
         % subgroups of faces that make cut-cells.

         % Make list of edges:
         EDGFAC_CELL = IBM_UNDEFINED*ones(SIZE_CFELEM_EDGFAC,SIZE_CEELEM_EDGFAC);
         FACEDG_CELL = IBM_UNDEFINED*ones(SIZE_CEELEM_FACEDG,SIZE_CFELEM_FACEDG);

         % Here reallocate FACEDG_CELL if NFACE_CELL > SIZE_CFELEM_FACEDG:
%          if (NFACE_CELL > SIZE_CFELEM_FACEDG)
%             DFCT = CEILING(REAL(NFACE_CELL-SIZE_CFELEM_FACEDG,EB)/REAL(DELTA_FACE,EB))
%             ALLOCATE(FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG+DFCT*DELTA_FACE));
%             FACEDG_CELL_AUX = IBM_UNDEFINED
%             % Copy data into FACEDG_CELL_AUX:
%             FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG) = &
%                 FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG)
%             % New SIZE_CFELEM_FACEDG:
%             SIZE_CFELEM_FACEDG = SIZE_CFELEM_FACEDG + DFCT*DELTA_FACE
%             DEALLOCATE(FACEDG_CELL); ALLOCATE(FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG))
%             % Dump data back into FACEDG_CELL:
%             FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG) = &
%             FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG)
%             DEALLOCATE(FACEDG_CELL_AUX)
%          end

         for IFACE=1:NFACE_CELL
            NIEDGE = FACE_CELL(1,IFACE);

            % Here reallocate if NIEDGE > SIZE_CEELEM_FACEDG:
%             if (NIEDGE > SIZE_CEELEM_FACEDG)
%                DFCT = CEILING(REAL(NIEDGE-SIZE_CEELEM_FACEDG,EB)/REAL(DELTA_EDGE,EB))
%                ALLOCATE(FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG+DFCT*DELTA_EDGE,1:SIZE_CFELEM_FACEDG));
%                FACEDG_CELL_AUX = IBM_UNDEFINED
%                % Copy data into FACEDG_CELL_AUX:
%                FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG) = &
%                    FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG)
%                % New SIZE_CEELEM_FACEDG:
%                SIZE_CEELEM_FACEDG = SIZE_CEELEM_FACEDG + DFCT*DELTA_EDGE
%                DEALLOCATE(FACEDG_CELL); ALLOCATE(FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG))
%                % Dump data back into FACEDG_CELL:
%                FACEDG_CELL(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG) = &
%                FACEDG_CELL_AUX(1:SIZE_CEELEM_FACEDG,1:SIZE_CFELEM_FACEDG)
%                DEALLOCATE(FACEDG_CELL_AUX)
%                DEALLOCATE(IPTS); ALLOCATE(IPTS(1:SIZE_CEELEM_FACEDG+1))
%             end

            IPTS(1:NIEDGE) = FACE_CELL(2:NIEDGE+1,IFACE)'; IPTS(NIEDGE+1) = FACE_CELL(2,IFACE);
            for IEDGE=1:NIEDGE
               SEG(NOD1:NOD2)= [ IPTS(IEDGE), IPTS(IEDGE+1) ];
               INLIST = false;
               for ISEG=1:NSEG_CELL
                  TEST1 = (SEG_CELL(NOD1,ISEG) == SEG(NOD1)) && (SEG_CELL(NOD2,ISEG) == SEG(NOD2));
                  TEST2 = (SEG_CELL(NOD2,ISEG) == SEG(NOD1)) && (SEG_CELL(NOD1,ISEG) == SEG(NOD2));

                  if ( TEST1 || TEST2 )
                     INLIST = true;
                     break
                  end
               end
               if (~INLIST)
                  NSEG_CELL = NSEG_CELL + 1;

                  % Test the NSEG_CELL doesn't overrun SIZE_CEELEM_EDGFAC, if so reallocate EDGFAC_CELL:
%                   if(NSEG_CELL > SIZE_CEELEM_EDGFAC)
%                      % 1. EDGFAC_CELL:
%                      ALLOCATE(EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC+DELTA_EDGE));
%                      EDGFAC_CELL_AUX = IBM_UNDEFINED
%                      % Copy data into EDGFAC_CELL_AUX:
%                      EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC) = &
%                          EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC)
%                      % 1. SEG_CELL:
%                      ALLOCATE(SEG_CELL_AUX(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC+DELTA_EDGE)); SEG_CELL_AUX = IBM_UNDEFINED
%                      % Copy data to SEG_CELL_AUX:
%                      SEG_CELL_AUX(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC) = SEG_CELL(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC)
% 
%                      % New SIZE_CEELEM_EDGFAC:
%                      SIZE_CEELEM_EDGFAC = SIZE_CEELEM_EDGFAC + DELTA_EDGE
% 
%                      % 2. EDGFAC_CELL:
%                      DEALLOCATE(EDGFAC_CELL); ALLOCATE(EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC))
%                      % Dump data back into EDGFAC_CELL:
%                      EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC) = &
%                      EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC)
%                      DEALLOCATE(EDGFAC_CELL_AUX)
%                      % 2. SEG_CELL:
%                      DEALLOCATE(SEG_CELL); ALLOCATE(SEG_CELL(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC))
%                      % Dump data back into SEG_CELL:
%                      SEG_CELL(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC) = SEG_CELL_AUX(NOD1:NOD2,1:SIZE_CEELEM_EDGFAC)
%                      DEALLOCATE(SEG_CELL_AUX)
%                   end
                  SEG_CELL(NOD1:NOD2,NSEG_CELL) = SEG(NOD1:NOD2)';
                  NEF = 1;
                  EDGFAC_CELL(1,NSEG_CELL)    =       NEF;
                  EDGFAC_CELL(NEF+1,NSEG_CELL)=     IFACE;
                  FACEDG_CELL(IEDGE,IFACE)    = NSEG_CELL;
               else
                  NEF = EDGFAC_CELL(1,ISEG) + 1;
                  % Test NEF+1 doesn't overrun SIZE_CFELEM_EDGFAC, if so reallocate EDGFAC_CELL:
%                   if(NEF+1 > SIZE_CFELEM_EDGFAC)
%                      ALLOCATE(EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC+DELTA_FACE,1:SIZE_CEELEM_EDGFAC));
%                      EDGFAC_CELL_AUX = IBM_UNDEFINED
%                      % Copy data into EDGFAC_CELL_AUX:
%                      EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC) = &
%                          EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC)
%                      % New SIZE_CFELEM_EDGFAC:
%                      SIZE_CFELEM_EDGFAC = SIZE_CFELEM_EDGFAC + DELTA_FACE
%                      DEALLOCATE(EDGFAC_CELL); ALLOCATE(EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC))
%                      % Dump data back into EDGFAC_CELL:
%                      EDGFAC_CELL(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC) = &
%                      EDGFAC_CELL_AUX(1:SIZE_CFELEM_EDGFAC,1:SIZE_CEELEM_EDGFAC)
%                      DEALLOCATE(EDGFAC_CELL_AUX)
%                   end
                  EDGFAC_CELL(1,ISEG)         =   NEF;
                  EDGFAC_CELL(NEF+1,ISEG)     = IFACE;
                  FACEDG_CELL(IEDGE,IFACE)    =  ISEG;
               end
            end
         end

         % Then  loop is on faces that have all regular edges,
         % that is, edges shared with only one another face:
         % Reallocate FACECELL_NUM if NFACE_CELL > SIZE(FACECELL_NUM,DIM=1):
%          NUM_FACE = SIZE(FACECELL_NUM,DIM=1)
%          if (NFACE_CELL > NUM_FACE)
%             DFCT = CEILING(REAL(NFACE_CELL-NUM_FACE,EB)/REAL(DELTA_FACE,EB))
%             DEALLOCATE(FACECELL_NUM); ALLOCATE(FACECELL_NUM(1:NFACE_CELL+DFCT*DELTA_FACE))
%          end

         FACECELL_NUM = zeros(1,NFACE_CELL);
         ICELL        = 1;
         IFACE        = 1;
         NUM_FACE     = NFACE_CELL;
         CTVAL2       = 0;
         MAXSEG       = max(FACE_CELL(1,1:NFACE_CELL));
         THRES        = (MAXSEG*NFACE_CELL)^2;
         CYCLE_CELL   =false;
         % Now double infinite loops:
         while(1)
            CTVAL        = 0;
            while(1)

               NEWFACE = false;
               NFACEI  = FACE_CELL(1,IFACE);

               % Now loop to find new face:
               for ISEG=1:NFACEI
                  LOCSEG = FACEDG_CELL(ISEG,IFACE);
                  if ( EDGFAC_CELL(1,LOCSEG) == 2 ) % Found a regular edge
                     for JJ=2:EDGFAC_CELL(1,LOCSEG)+1
                        JFACE = EDGFAC_CELL(JJ,LOCSEG);
                        % Drop for same face:
                        if ( IFACE == JFACE ); continue; end
                        % Drop if face already counted:
                        if ( FACECELL_NUM(JFACE) > 0 ); continue; end

                        % New face, not counted:
                        FACECELL_NUM(JFACE) =  ICELL;
                        NEWFACE             = true;
                        NUM_FACE            = NUM_FACE-1;
                        break
                     end
                  end
                  if (NEWFACE)
                     IFACE = JFACE;
                     break
                  end
               end

               % Test for all faces that have regular edges with faces that belong to icell:
               if (~NEWFACE)
                  for KFACE=1:NFACE_CELL
                     KFACEFLG=false;
                     if ( FACECELL_NUM(KFACE) == 0 ) % Not associated yet
                        NFACEK = FACE_CELL(1,KFACE);
                        for ISEG=1:NFACEK
                           LOCSEG = FACEDG_CELL(ISEG,KFACE);
                           if ( EDGFAC_CELL(1,LOCSEG) == 2) % Found a regular edge
                              for JJ=2:EDGFAC_CELL(1,LOCSEG)+1
                                 JFACE = EDGFAC_CELL(JJ,LOCSEG);
                                 if ( KFACE == JFACE ); continue; end
                                 if ( FACECELL_NUM(JFACE) ~= ICELL); continue; end
                                 % New face, not counted:
                                 FACECELL_NUM(KFACE) = FACECELL_NUM(JFACE);
                                 NEWFACE             = true;
                                 IFACE               = KFACE;
                                 NUM_FACE            = NUM_FACE-1;
                                 KFACEFLG            = true;
                                 break
                              end
                           end
                           if(KFACEFLG); break; end
                        end
                     end
                     if(KFACEFLG); break; end
                  end
               end

               % Haven't found new face, either num_face=0, or we need a new icell:
               if (~NEWFACE); break; end
               
               CTVAL = CTVAL + 1;
               if CTVAL > THRES
                   disp(['Inner Special cell = ' num2str(I) ',' num2str(J) ',' num2str(K)])
                   CYCLE_CELL=true;
                   break
               end

            end
            % Test if there are any faces left:
            if ( NUM_FACE <= 0 )
               break
            else % New cell, find new face set iface
               FNDFLG = 0;
               for IFACE=1:NFACE_CELL
                  if (FACECELL_NUM(IFACE) == 0) % NOT COUNTED YET.
                       % ASSUMES IT HAS AT LEAST ONE REGULAR EDGE.
                       ICELL = ICELL + 1;
                       FNDFLG=1;
                       break
                   end
               end
               if(~FNDFLG); break; end % Case all faces associated.
            end
            CTVAL2 = CTVAL2 + 1;
            if CTVAL2 > THRES
                disp(['Outer Special cell = ' num2str(I) ',' num2str(J) ',' num2str(K) ...
                      ', NFACE_CELL=' num2str(NFACE_CELL) ', NUM_FACE=' num2str(NUM_FACE)])
                CYCLE_CELL=true;
            end
            if(CYCLE_CELL); break; end
         end

         if(CYCLE_CELL)
             CELLRT(I,J,K) = true;
             MESHES(NM).N_SPCELL = MESHES(NM).N_SPCELL+1;
             MESHES(NM).SPCELL_LIST(IAXIS:KAXIS,MESHES(NM).N_SPCELL) = [I J K]';
             
             IDCF = MESHES(NM).CCVAR(I,J,K,IBM_IDCF);
             disp(['I,J,K,IDCF=' num2str([I J K IDCF])])
             NIBFACE    = 0;
             NFACE_CELL = N_GAS_CFACES + NIBFACE;
             if (IDCF > 0)
                 
                 IBOD = 1; ITRI = 1;
                 if (MESHES(NM).CUT_FACE(IDCF).NFACE > 0)
                   IBOD = MESHES(NM).CUT_FACE(IDCF).BODTRI(1,1);
                   ITRI = MESHES(NM).CUT_FACE(IDCF).BODTRI(2,1);
                 end
                 
                 NIBFACE    = 0;
                 XYZVERT    = zeros(KAXIS,IBM_MAXVERTS_CELL);
                 NVERT_CELL = 0;
                 CFELEM     = zeros(5,1);
                 % Define from SOLID FACES CFACES for the cell:
                 IED = I-FCELL; JED = J-FCELL; KED = K-FCELL;
                 for MYAXIS=IAXIS:KAXIS
                     switch(MYAXIS)
                         case(IAXIS)
                             XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD1,LOW_IND)  = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,LOW_IND)  = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,LOW_IND)  = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,LOW_IND)  = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';
                             AREAI    = DYCELL(J) * DZCELL(K);
                         case(JAXIS)
                             XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD1,LOW_IND)  = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,LOW_IND)  = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,LOW_IND)  = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,LOW_IND)  = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';
                             AREAI    = DXCELL(I) * DZCELL(K);
                         case(KAXIS)
                             XYZLH(IAXIS:KAXIS,NOD1,HIGH_IND) = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,HIGH_IND) = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,HIGH_IND) = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,HIGH_IND) = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED+1) ]';
                             XYZLH(IAXIS:KAXIS,NOD1,LOW_IND)  = [ XFACE(IED  ), YFACE(JED  ), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD2,LOW_IND)  = [ XFACE(IED+1), YFACE(JED  ), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD3,LOW_IND)  = [ XFACE(IED+1), YFACE(JED+1), ZFACE(KED  ) ]';
                             XYZLH(IAXIS:KAXIS,NOD4,LOW_IND)  = [ XFACE(IED  ), YFACE(JED+1), ZFACE(KED  ) ]';
                             AREAI    = DXCELL(I) * DYCELL(J);
                     end                     
                     
                     for SIDE=LOW_IND:HIGH_IND
                        if (FSID_XYZ(SIDE ,MYAXIS) ~= IBM_SOLID); continue; end
                        NIBFACE = NIBFACE + 1;                       
                        
                        % Define vertices of CFACE and insert add to
                        % MESHES(NM).CUT_FACE(IDCF).XYZVERT
                        NP = 0; 
                        XYZC(IAXIS:KAXIS) = 0.;
                        for IP=NOD1:NOD4
                           % xl,yl,zl
                           XYZ(IAXIS:KAXIS) = XYZLH(IAXIS:KAXIS,IP,SIDE)';
                           XYZC(IAXIS:KAXIS)= XYZC(IAXIS:KAXIS) + XYZ(IAXIS:KAXIS);
                           [NVERT_CELL,INOD,XYZVERT]=INSERT_FACE_VERT_LOC(IBM_MAXVERTS_CELL,XYZ,NVERT_CELL,XYZVERT);
                           NP = NP + 1;
                           CFELEM(1)    = NP;
                           CFELEM(NP+1) = INOD;
                        end               
                        
                        % Define CFELEM connectivity, also CFACE area and
                        % Centroid add to corresponding CUT_FACE(IDCF)
                        % entries.
                        MESHES(NM).CUT_FACE(IDCF).CFELEM(1:5,NIBFACE) = CFELEM(1:5)';
                        MESHES(NM).CUT_FACE(IDCF).AREA(NIBFACE) = AREAI;
                        MESHES(NM).CUT_FACE(IDCF).XYZCEN(IAXIS:KAXIS,NIBFACE) = 0.25*XYZC(IAXIS:KAXIS)';
                        % Fields for cut-cell volume/centroid computation:
                        % dot(i,nc)*int(x)dA:
                        MESHES(NM).CUT_FACE(IDCF).INXAREA(NIBFACE)   = 0.;
                        % dot(i,nc)*int(x^2)dA:
                        MESHES(NM).CUT_FACE(IDCF).INXSQAREA(NIBFACE) = 0.;
                        % dot(j,nc)*int(y^2)dA:
                        MESHES(NM).CUT_FACE(IDCF).JNYSQAREA(NIBFACE) = 0.;
                        % dot(k,nc)*int(z^2)dA:
                        MESHES(NM).CUT_FACE(IDCF).KNZSQAREA(NIBFACE) = 0.;

                        % Define Body-triangle reference: 
                        MESHES(NM).CUT_FACE(IDCF).BODTRI(1:2,NIBFACE)= [ IBOD, ITRI ];

                        % Assign surf-index: Depending on GEOMETRY:
                        % Here we might just add the INERT SURF_ID:
                        MESHES(NM).CUT_FACE(IDCF).SURF_INDEX(NIBFACE) = GEOM(IBOD).SURFS(ITRI);
                        
                        
                        
                        
                        % Finally add to FACE_LIST from N_GAS_CFACES on:
                        NFACE_CELL = N_GAS_CFACES + NIBFACE;
                        FACE_LIST(1:IBM_NPARAM_CCFACE,NFACE_CELL) = [ IBM_FTYPE_CFINB, 0, 0, IDCF, NIBFACE ]';                    
                        
                     
                     end
                 end
                 MESHES(NM).CUT_FACE(IDCF).NFACE = NIBFACE;             
                 MESHES(NM).CUT_FACE(IDCF).XYZVERT(IAXIS:KAXIS,1:NVERT_CELL) = XYZVERT(IAXIS:KAXIS,1:NVERT_CELL);
             end
             
             % Now define a coarse cut-cell (no INBOUNDARY cut-faces):
             NCELL      = 1;
             %NFACE_CELL = N_GAS_CFACES+NIBFACE;
             CCELEM(1:NFACE_CELL+1) = [NFACE_CELL; [1:NFACE_CELL]'];
             VOL = zeros(1,NCELL);
             XYZCEN = zeros(KAXIS,NCELL);
             VOL(NCELL) = DXCELL(I)*DYCELL(J)*DZCELL(K);
             XYZCEN(IAXIS:KAXIS,NCELL) = [XCELL(I) YCELL(J) ZCELL(K)]';


         else

         % Create CCELEM array:
         NCELL = max(FACECELL_NUM(:));
%          % Test NCELL not > SIZE_CELL_CCELEM; NFACE_CELL not > SIZE_FACE_CCELEM:
%          if (NFACE_CELL > SIZE_FACE_CCELEM)
%             DFCT = CEILING(REAL(NFACE_CELL-SIZE_FACE_CCELEM,EB)/REAL(DELTA_FACE,EB))
%             SIZE_FACE_CCELEM = SIZE_FACE_CCELEM + DFCT*DELTA_FACE
%             DEALLOCATE(CCELEM)
%             ALLOCATE(CCELEM(1:SIZE_FACE_CCELEM+1,1:SIZE_CELL_CCELEM))
%          end
%          if (NCELL > SIZE_CELL_CCELEM)
%             DFCT = CEILING(REAL(NCELL-SIZE_CELL_CCELEM,EB)/REAL(DELTA_CELL,EB))
%             SIZE_CELL_CCELEM = SIZE_CELL_CCELEM + DFCT*DELTA_CELL
%             DEALLOCATE(CCELEM,VOL,XYZCEN)
%             ALLOCATE(CCELEM(1:SIZE_FACE_CCELEM+1,1:SIZE_CELL_CCELEM))
%             ALLOCATE(VOL(1:SIZE_CELL_CCELEM),XYZCEN(IAXIS:KAXIS,1:SIZE_CELL_CCELEM))
%          end
         CCELEM= IBM_UNDEFINED*ones(SIZE_FACE_CCELEM+1,SIZE_CELL_CCELEM);
         for ICELL=1:NCELL
            NP = 0;
            for IFACE=1:NFACE_CELL
               if ( FACECELL_NUM(IFACE) == ICELL )
                  NP = NP + 1;
                  CCELEM(1,ICELL)    =    NP;
                  CCELEM(NP+1,ICELL) = IFACE;
               end
            end
         end

         % Compute volumes and centroids for the found cut-cells:
         VOL = zeros(1,NCELL);
         XYZCEN = zeros(KAXIS,NCELL);
         for ICELL=1:NCELL
            NP = CCELEM(1,ICELL);
            for II=2:NP+1
               IFACE = CCELEM(II,ICELL);

               % Volume:
               VOL(ICELL) = VOL(ICELL) + AREAVARS(1,IFACE);

               % xyzcen:
               XYZCEN(IAXIS:KAXIS,ICELL) = XYZCEN(IAXIS:KAXIS,ICELL)+AREAVARS(2:4,IFACE);
            end
            if(VOL(ICELL) < GEOMEPS) % Volume too small for correct calculation of XYZCEN-> take cartcell centroid.
               JJ = 0;
               VOL(ICELL) = abs(VOL(ICELL));
               XYZCEN(IAXIS:KAXIS,ICELL) = [ XCELL(I), YCELL(J), ZCELL(K) ]';
            else
               % divide xyzcen by 2*vol:
               XYZCEN(IAXIS:KAXIS,ICELL) = XYZCEN(IAXIS:KAXIS,ICELL) / (2.*VOL(ICELL));
            end
         end

         end

         % Load into CUT_CELL data structure
         NCUTCELL = MESHES(NM).N_CUTCELL_MESH + MESHES(NM).N_GCCUTCELL_MESH + 1;
         if (IBNDINT==LOW_IND)
            MESHES(NM).N_CUTCELL_MESH   = NCUTCELL;
         else
            MESHES(NM).N_GCCUTCELL_MESH = MESHES(NM).N_GCCUTCELL_MESH + 1;
         end
         MESHES(NM).CCVAR(I,J,K,IBM_IDCC)            = NCUTCELL;

         % Resize array MESHES(NM).CUT_CELL if necessary:
         % CALL CUT_CELL_ARRAY_REALLOC(NM,NCUTCELL)

         % Add cut-cell NCUTCELL entry:
         MESHES(NM).CUT_CELL(NCUTCELL).IJK(IAXIS:KAXIS) = [ I, J, K ];
         MESHES(NM).CUT_CELL(NCUTCELL).NCELL     = NCELL;
         MESHES(NM).CUT_CELL(NCUTCELL).NFACE_CELL= NFACE_CELL;
         NCFACE_CUTCELL = max(CCELEM(1,1:NCELL)) + 1;
         %CALL NEW_CELL_ALLOC(NM,NCUTCELL,NCELL,NFACE_CELL,NCFACE_CUTCELL)         
         MESHES(NM).CUT_CELL(NCUTCELL).CCELEM(1:NCFACE_CUTCELL,1:NCELL) = CCELEM(1:NCFACE_CUTCELL,1:NCELL);
         MESHES(NM).CUT_CELL(NCUTCELL).FACE_LIST(1:IBM_NPARAM_CCFACE,1:NFACE_CELL) = ...
         FACE_LIST(1:IBM_NPARAM_CCFACE,1:NFACE_CELL);
         MESHES(NM).CUT_CELL(NCUTCELL).VOLUME(1:NCELL)                  = VOL(1:NCELL);
         MESHES(NM).CUT_CELL(NCUTCELL).XYZCEN(IAXIS:KAXIS,1:NCELL)      = XYZCEN(IAXIS:KAXIS,1:NCELL);

      end % I
   end % J
end % K

end

ierr=0;

return