close all 
clear all
clc

% Declarations, call set constants:
global MAX_DIM IAXIS JAXIS KAXIS NOD1 NOD2 NOD3 NOD4 
global NMESHES MESHES BODINT_PLANE
global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_BDNUM_CRS_AUX IBM_IS_CRS2_AUX
global NM N_GEOMETRY GEOM IBM_INBOUNDARY IBM_INBOUNDCC IBM_INBOUNDCF IBM_CGSC IBM_CUTCFE
global IBM_IDCF IBM_IDCE IBM_FGSC IBM_IDCC IBM_GASPHASE
global XCELL DXCELL XFACE DXFACE 
global YCELL DYCELL YFACE DYFACE
global ZCELL DZCELL ZFACE DZFACE
global IBM_FTYPE_RGGAS IBM_FTYPE_CFGAS IBM_FTYPE_CFINB

global LOW_IND HIGH_IND

global basedir CELLRT

global BODINT_PLANE3 GEOMEPS

global XAXIS XRAY

global SEG_FACE_F SEG_FACE2_F XYZVERT_F NFACE_F CFELEM_F NVERT_F

plot_cutedges=false;

[ierr]=SET_CONSTANTS();

% Case name and directory:
addpath ../
basedir = '/Users/mnv/Documents/FIREMODELS_FORK/fds/GEOM_Intersection/';
%casename='2A_V1_cat';
%casename='geom_self_intersection';
%casename='leak_test_4';
%casename= 'geom_rad_2'; 
%casename='geom_intersect'; 
%casename='two_spheres';
casename='geom_extruded_poly';

% Load and plot Geometries:
[GEOM,N_GEOMETRY]=load_geometries(basedir,casename);
newfig =1;
[ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig);


% Load and plot Meshes:
[MESHES,NMESHES]=load_meshes(basedir,casename);
newfig=0;
[ierr]=plot_meshes(NMESHES,MESHES,newfig);

tstart = tic;
[ierr]=SET_CUTCELLS_3D(basedir,casename,plot_cutedges);
toc(tstart)


figure
hold on
axis equal; box on;
view([45 45])
xlabel('X')
ylabel('Y')
zlabel('Z')
 NM=1;
a = MESHES(NM).DX(6)/2;
% Plot Boundary cut-faces:
for ICF=1:MESHES(NM).N_CUTFACE_MESH
   if(MESHES(NM).CUT_FACE(ICF).STATUS ~= IBM_INBOUNDARY); continue; end
   NFACE=MESHES(NM).CUT_FACE(ICF).NFACE;
   XYZVERT = MESHES(NM).CUT_FACE(ICF).XYZVERT;
   for JCF=1:NFACE
       NELEM  = MESHES(NM).CUT_FACE(ICF).CFELEM(1,JCF);
       CFELEM = MESHES(NM).CUT_FACE(ICF).CFELEM(2:NELEM+1,JCF);
       
       IG     = MESHES(NM).CUT_FACE(ICF).BODTRI(1,JCF);
       
       if(IG==1)
       [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'r');
       else
       [hp]=patch(XYZVERT(IAXIS,CFELEM),XYZVERT(JAXIS,CFELEM),XYZVERT(KAXIS,CFELEM),'b');
       end
       set(hp,'FaceAlpha',0.3) %,'EdgeAlpha',0.4) % .4  
       
%        X0 = 1/NELEM*sum(XYZVERT(IAXIS,CFELEM));
%        Y0 = 1/NELEM*sum(XYZVERT(JAXIS,CFELEM));
%        Z0 = 1/NELEM*sum(XYZVERT(KAXIS,CFELEM));
%        TRI= MESHES(NM).CUT_FACE(ICF).BODTRI(2,JCF);     
%        N  = GEOM(IG).FACES_NORMAL(IAXIS:KAXIS,TRI)';
%        plot3([X0 X0+a*N(IAXIS)],[Y0 Y0+a*N(JAXIS)],[Z0 Z0+a*N(KAXIS)],'k')
       
   end
   
%    I = MESHES(NM).CUT_FACE(ICF).IJK(IAXIS);
%    J = MESHES(NM).CUT_FACE(ICF).IJK(JAXIS);
%    K = MESHES(NM).CUT_FACE(ICF).IJK(KAXIS);
%    
%    P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
%    P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
%    P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
%    P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
%    P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
%    P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
%    P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
%    P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];
%    
%    FC(1,:) = [ 1 4 8 5 ]; VR(1)=MESHES(NM).FCVAR(I-1,J  ,K  ,IBM_FGSC,IAXIS);
%    FC(2,:) = [ 2 3 7 6 ]; VR(2)=MESHES(NM).FCVAR(I  ,J  ,K  ,IBM_FGSC,IAXIS);
%    FC(3,:) = [ 1 5 6 2 ]; VR(3)=MESHES(NM).FCVAR(I  ,J-1,K  ,IBM_FGSC,JAXIS);
%    FC(4,:) = [ 4 8 7 3 ]; VR(4)=MESHES(NM).FCVAR(I  ,J  ,K  ,IBM_FGSC,JAXIS);
%    FC(5,:) = [ 1 2 3 4 ]; VR(5)=MESHES(NM).FCVAR(I  ,J  ,K-1,IBM_FGSC,KAXIS);
%    FC(6,:) = [ 5 6 7 8 ]; VR(6)=MESHES(NM).FCVAR(I  ,J  ,K  ,IBM_FGSC,KAXIS);
%    
%    for IFC=1:6
%        if(VR(IFC)==IBM_GASPHASE); continue; end
%        [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'m');
%        set(hp,'FaceAlpha',0.2)
%    end


end


% % Plot special cells:
disp(['N_SPCELL=' num2str(MESHES(NM).N_SPCELL)])
% for ICELL=1:MESHES(NM).N_SPCELL
%     I=MESHES(NM).SPCELL_LIST(IAXIS,ICELL);
%     J=MESHES(NM).SPCELL_LIST(JAXIS,ICELL);
%     K=MESHES(NM).SPCELL_LIST(KAXIS,ICELL);
%     
%     P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
%     P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
%     P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
%     P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
%     P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
%     P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
%     P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
%     P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];
% 
%     FC(1,:) = [ 1 4 8 5 ];
%     FC(2,:) = [ 2 3 7 6 ];
%     FC(3,:) = [ 1 5 6 2 ];
%     FC(4,:) = [ 4 8 7 3 ];
%     FC(5,:) = [ 1 2 3 4 ];
%     FC(6,:) = [ 5 6 7 8 ];
% 
%     for IFC=1:6
%     [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'r');
%        set(hp,'FaceAlpha',0.3)
%     end
% end

% for K=MESHES(NM).KLO_CELL:MESHES(NM).KHI_CELL
%    for J=MESHES(NM).JLO_CELL:MESHES(NM).JHI_CELL
%       for I=MESHES(NM).ILO_CELL:MESHES(NM).IHI_CELL
% 
%         if(MESHES(NM).CCVAR(I,J,K,IBM_CGSC) ~= IBM_CUTCFE); continue; end
%         if(~CELLRT(I,J,K)); continue; end
%         P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
%         P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
%         P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
%         P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
%         P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
%         P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
%         P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
%         P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];
%         
%         FC(1,:) = [ 1 4 8 5 ];
%         FC(2,:) = [ 2 3 7 6 ];
%         FC(3,:) = [ 1 5 6 2 ];
%         FC(4,:) = [ 4 8 7 3 ];
%         FC(5,:) = [ 1 2 3 4 ];
%         FC(6,:) = [ 5 6 7 8 ];
%         
%         for IFC=1:6
%             [hp]=patch(P(FC(IFC,:),IAXIS),P(FC(IFC,:),JAXIS),P(FC(IFC,:),KAXIS),'b');
%             set(hp,'FaceAlpha',0.3)
%         end
%   
%       end
%    end
% end

return

figure
hold on
axis equal; box on;
view([45 45])
xlabel('X')
ylabel('Y')
zlabel('Z')

% for ICC2=1:MESHES(NM).N_SPCELL    
%    ICC = MESHES(NM).CCVAR(MESHES(NM).SPCELL_LIST(IAXIS,ICC2),...
%                           MESHES(NM).SPCELL_LIST(JAXIS,ICC2),...
%                           MESHES(NM).SPCELL_LIST(KAXIS,ICC2),IBM_IDCC);
  
for ICC=1:MESHES(NM).N_CUTCELL_MESH+MESHES(NM).N_GCCUTCELL_MESH
   I      = MESHES(NM).CUT_CELL(ICC).IJK(IAXIS);
   J      = MESHES(NM).CUT_CELL(ICC).IJK(JAXIS);
   K      = MESHES(NM).CUT_CELL(ICC).IJK(KAXIS);
   NCELL  = MESHES(NM).CUT_CELL(ICC).NCELL;
   
   P(1,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K-1)];
   P(2,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K-1)];
   P(3,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K-1)];
   P(4,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K-1)];
   P(5,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J-1) ZFACE(K  )];
   P(6,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J-1) ZFACE(K  )];
   P(7,IAXIS:KAXIS) = [XFACE(I  ) YFACE(J  ) ZFACE(K  )];
   P(8,IAXIS:KAXIS) = [XFACE(I-1) YFACE(J  ) ZFACE(K  )];

   for JCC=1:NCELL
      NFCELL = MESHES(NM).CUT_CELL(ICC).CCELEM(1,JCC);
      
      for ICCF=1:NFCELL
          IFACE=MESHES(NM).CUT_CELL(ICC).CCELEM(ICCF+1,JCC);
          
          switch(MESHES(NM).CUT_CELL(ICC).FACE_LIST(1,IFACE))
              case(IBM_FTYPE_RGGAS) % REGULAR GASPHASE
              LOWHIGH = MESHES(NM).CUT_CELL(ICC).FACE_LIST(2,IFACE);
              X1AXIS  = MESHES(NM).CUT_CELL(ICC).FACE_LIST(3,IFACE);
              ILH     = LOWHIGH - 1;
              switch(X1AXIS)
                  case(IAXIS)
                      if(LOWHIGH== LOW_IND); FC = [ 1 4 8 5 ]; end
                      if(LOWHIGH==HIGH_IND); FC = [ 2 3 7 6 ]; end
                  case(JAXIS)
                      if(LOWHIGH== LOW_IND); FC = [ 1 5 6 2 ]; end
                      if(LOWHIGH==HIGH_IND); FC = [ 4 8 7 3 ]; end
                  case(KAXIS)
                      if(LOWHIGH== LOW_IND); FC = [ 1 2 3 4 ]; end
                      if(LOWHIGH==HIGH_IND); FC = [ 5 6 7 8 ]; end
              end

              [hp]=patch(P(FC(:),IAXIS),P(FC(:),JAXIS),P(FC(:),KAXIS),'b');
              set(hp,'FaceAlpha',0.2,'EdgeAlpha',0.3)

              case(IBM_FTYPE_CFGAS)                  
              IFC2    = MESHES(NM).CUT_CELL(ICC).FACE_LIST(4,IFACE);
              IFACE2  = MESHES(NM).CUT_CELL(ICC).FACE_LIST(5,IFACE);

              PP      = MESHES(NM).CUT_FACE(IFC2).XYZVERT';
              NP      = MESHES(NM).CUT_FACE(IFC2).CFELEM(1,IFACE2);
              FC      = MESHES(NM).CUT_FACE(IFC2).CFELEM(2:NP+1,IFACE2);
              
              [hp]=patch(PP(FC(:),IAXIS),PP(FC(:),JAXIS),PP(FC(:),KAXIS),'b');
              set(hp,'FaceAlpha',0.2,'EdgeAlpha',0.3)

              case(IBM_FTYPE_CFINB)
              IFC2    = MESHES(NM).CUT_CELL(ICC).FACE_LIST(4,IFACE);
              IFACE2  = MESHES(NM).CUT_CELL(ICC).FACE_LIST(5,IFACE);

              PP      = MESHES(NM).CUT_FACE(IFC2).XYZVERT';
              NP      = MESHES(NM).CUT_FACE(IFC2).CFELEM(1,IFACE2);
              FC      = MESHES(NM).CUT_FACE(IFC2).CFELEM(2:NP+1,IFACE2);
              
              IG      = MESHES(NM).CUT_FACE(IFC2).BODTRI(1,IFACE2);
              
              if(IG==1)
              [hp]=patch(PP(FC(:),IAXIS),PP(FC(:),JAXIS),PP(FC(:),KAXIS),'r');                  
              else
              [hp]=patch(PP(FC(:),IAXIS),PP(FC(:),JAXIS),PP(FC(:),KAXIS),'k');
              end
              set(hp,'FaceAlpha',0.4,'EdgeAlpha',0.4)
              
%               % PLot Normal:
%               if (NP==3)
%               PP  = PP'; 
%               D12 = PP(IAXIS:KAXIS,FC(2))-PP(IAXIS:KAXIS,FC(1));
%               D23 = PP(IAXIS:KAXIS,FC(3))-PP(IAXIS:KAXIS,FC(2));
%               X0 = 1/NP*sum(PP(IAXIS,FC));
%               Y0 = 1/NP*sum(PP(JAXIS,FC));
%               Z0 = 1/NP*sum(PP(KAXIS,FC));
%             
%               N  = cross(D12,D23); N=N/norm(N);
%               plot3([X0 X0+a*N(IAXIS)],[Y0 Y0+a*N(JAXIS)],[Z0 Z0+a*N(KAXIS)],'k')
%               end

          end  
          
      end
        
   end
end

% Sanity test
for ICC=1:MESHES(NM).N_CUTCELL_MESH+MESHES(NM).N_GCCUTCELL_MESH
   I      = MESHES(NM).CUT_CELL(ICC).IJK(IAXIS);
   J      = MESHES(NM).CUT_CELL(ICC).IJK(JAXIS);
   K      = MESHES(NM).CUT_CELL(ICC).IJK(KAXIS);
   NCELL  = MESHES(NM).CUT_CELL(ICC).NCELL;
   for JCC=1:NCELL
       vol = MESHES(NM).CUT_CELL(ICC).VOLUME(JCC);
       if (MESHES(NM).CUT_CELL(ICC).VOLUME(JCC) < 0)
           disp(['Negative Volume: ' num2str(ICC) ', ' num2str(JCC) ', vol: ' num2str(vol)])
       end
       XYZCEN = MESHES(NM).CUT_CELL(ICC).XYZCEN(IAXIS:KAXIS,JCC)';
       if ( (XYZCEN(IAXIS) < XFACE(I-1)) || (XYZCEN(IAXIS) > XFACE(I)  ) )
           disp(['IAXIS XCEN outside: ' num2str(ICC) ', ' num2str(JCC) ', ' ...
                num2str([XFACE(I-1) XFACE(I) XYZCEN(IAXIS)]) ', CELLRT:' num2str(CELLRT(I,J,K))])
       end
       if ( (XYZCEN(JAXIS) < YFACE(J-1)) || (XYZCEN(JAXIS) > YFACE(J)  ) )
           disp(['JAXIS XCEN outside: ' num2str(ICC) ', ' num2str(JCC) ', ' ...
                num2str([YFACE(J-1) YFACE(J) XYZCEN(JAXIS)]) ', CELLRT:' num2str(CELLRT(I,J,K))])
       end
       if ( (XYZCEN(KAXIS) < ZFACE(K-1)) || (XYZCEN(KAXIS) > ZFACE(K)  ) )
           disp(['KAXIS XCEN outside: ' num2str(ICC) ', ' num2str(JCC) ', ' ...
                num2str([ZFACE(K-1) ZFACE(K) XYZCEN(KAXIS)]) ', CELLRT:' num2str(CELLRT(I,J,K))])
       end
   end
end


return
