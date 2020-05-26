function [MESHES,NMESHES]=load_meshes(basedir,casename)

global MAX_DIM NGUARD IAXIS KAXIS
global IBM_GASPHASE
global IBM_NVVARS IBM_VGSC IBM_NEVARS IBM_EGSC IBM_NFVARS IBM_FGSC IBM_NCVARS IBM_CGSC
global TRANS

[fid]=fopen([basedir casename '_meshes.dat'],'r');

NMESHES= str2num(fgetl(fid));

for NM=1:NMESHES
   DATA = str2num(fgetl(fid));
   MESHES(NM).IBAR = DATA( 2); 
   MESHES(NM).JBAR = DATA( 3); 
   MESHES(NM).KBAR = DATA( 4);
   MESHES(NM).XS   = DATA( 5);
   MESHES(NM).XF   = DATA( 6);
   MESHES(NM).YS   = DATA( 7);
   MESHES(NM).YF   = DATA( 8);
   MESHES(NM).ZS   = DATA( 9);
   MESHES(NM).ZF   = DATA(10);
   
   % X coordinate arrays:
   MESHES(NM).ISTR     = 1;
   MESHES(NM).ILO_FACE = NGUARD+0;
   MESHES(NM).ILO_CELL = NGUARD+1;
   MESHES(NM).IHI_FACE = MESHES(NM).ILO_FACE + MESHES(NM).IBAR;
   MESHES(NM).IHI_CELL = MESHES(NM).ILO_CELL + MESHES(NM).IBAR - 1;
   MESHES(NM).IEND     = MESHES(NM).IHI_FACE + NGUARD;   
   for I=MESHES(NM).ILO_FACE:MESHES(NM).IHI_FACE
      DATA = str2num(fgetl(fid));
      MESHES(NM).X(I)  = DATA(1);
      MESHES(NM).XC(I) = DATA(2);
      MESHES(NM).DXN(I)= DATA(3);
      MESHES(NM).DX(I) = DATA(4);
   end
   I = MESHES(NM).IHI_CELL + 1;
   MESHES(NM).XC(I) = MESHES(NM).XC(I-1) +  MESHES(NM).DX(I-1);
   MESHES(NM).DX(I) = MESHES(NM).DX(I-1);
   
   
   % Y Coordinate arrays:
   MESHES(NM).JSTR     = 1;
   MESHES(NM).JLO_FACE = NGUARD+0;
   MESHES(NM).JLO_CELL = NGUARD+1;
   MESHES(NM).JHI_FACE = MESHES(NM).JLO_FACE + MESHES(NM).JBAR;
   MESHES(NM).JHI_CELL = MESHES(NM).JLO_CELL + MESHES(NM).JBAR - 1;
   MESHES(NM).JEND     = MESHES(NM).JHI_FACE + NGUARD;   
   for J=MESHES(NM).JLO_FACE:MESHES(NM).JHI_FACE
      DATA = str2num(fgetl(fid));
      MESHES(NM).Y(J)  = DATA(1);
      MESHES(NM).YC(J) = DATA(2);
      MESHES(NM).DYN(J)= DATA(3);
      MESHES(NM).DY(J) = DATA(4);
   end   
   J = MESHES(NM).JHI_CELL + 1;
   MESHES(NM).YC(J) = MESHES(NM).YC(J-1) +  MESHES(NM).DY(J-1);
   MESHES(NM).DY(J) = MESHES(NM).DY(J-1);
   
   % Z Coordinate arrays:
   MESHES(NM).KSTR     = 1;
   MESHES(NM).KLO_FACE = NGUARD+0;
   MESHES(NM).KLO_CELL = NGUARD+1;
   MESHES(NM).KHI_FACE = MESHES(NM).KLO_FACE + MESHES(NM).KBAR;
   MESHES(NM).KHI_CELL = MESHES(NM).KLO_CELL + MESHES(NM).KBAR - 1;
   MESHES(NM).KEND     = MESHES(NM).KHI_FACE + NGUARD;   
   for K=MESHES(NM).KLO_FACE:MESHES(NM).KHI_FACE
      DATA = str2num(fgetl(fid));
      MESHES(NM).Z(K)  = DATA(1);
      MESHES(NM).ZC(K) = DATA(2);
      MESHES(NM).DZN(K)= DATA(3);
      MESHES(NM).DZ(K) = DATA(4);
   end   
   K = MESHES(NM).KHI_CELL + 1;
   MESHES(NM).ZC(K) = MESHES(NM).ZC(K-1) +  MESHES(NM).DZ(K-1);
   MESHES(NM).DZ(K) = MESHES(NM).DZ(K-1);
   
   
   
   
   % Define CC Mesh arrays:
   IEND = MESHES(NM).IEND;
   JEND = MESHES(NM).JEND;
   KEND = MESHES(NM).KEND;
   
   MESHES(NM).VERTVAR = zeros(IEND, JEND, KEND, IBM_NVVARS);
   MESHES(NM).VERTVAR(:,:,:,IBM_VGSC) = IBM_GASPHASE;
   
   MESHES(NM).ECVAR = zeros(IEND, JEND, KEND, IBM_NEVARS, MAX_DIM);
   MESHES(NM).ECVAR(:,:,:,IBM_EGSC,:) = IBM_GASPHASE;
   
   MESHES(NM).FCVAR = zeros(IEND, JEND, KEND, IBM_NFVARS, MAX_DIM);
   MESHES(NM).FCVAR(:,:,:,IBM_FGSC,:) = IBM_GASPHASE;
   
   MESHES(NM).CCVAR = zeros(IEND, JEND, KEND, IBM_NCVARS);
   MESHES(NM).CCVAR(:,:,:,IBM_CGSC) = IBM_GASPHASE;
   
   MESHES(NM).N_EDGE_CROSS   = 0;
   MESHES(NM).N_CUTEDGE_MESH = 0;
   MESHES(NM).N_CUTFACE_MESH = 0;
   MESHES(NM).N_CUTCELL_MESH = 0;
   
   MESHES(NM).N_GCCUTFACE_MESH = 0;
   MESHES(NM).N_GCCUTCELL_MESH = 0;
   
   TRANS(NM).NOC(IAXIS:KAXIS)= 0;
   
   
   
   DZMEAN=MESHES(NM).ZF-MESHES(NM).ZS/MESHES(NM).KBAR;
   if ( any( abs(MESHES(NM).DZ(MESHES(NM).KLO_CELL:MESHES(NM).KHI_CELL)-DZMEAN) > 1e-5*DZMEAN ) )
       TRANS(NM).NOC(KAXIS)= 1;
   end
   
   MESHES(NM).N_SPCELL = 0;
end

fclose(fid);

% DATA=load([basedir casename '_meshes.dat']);
% NMESHES=length(DATA(:,1));
% for NM=1:NMESHES
%    MESHES(NM).IBAR = DATA(NM, 2); 
%    MESHES(NM).JBAR = DATA(NM, 3); 
%    MESHES(NM).KBAR = DATA(NM, 4);
%    MESHES(NM).XS   = DATA(NM, 5);
%    MESHES(NM).XF   = DATA(NM, 6);
%    MESHES(NM).YS   = DATA(NM, 7);
%    MESHES(NM).YF   = DATA(NM, 8);
%    MESHES(NM).ZS   = DATA(NM, 9);
%    MESHES(NM).ZF   = DATA(NM,10);
%    MESHES(NM).DXYZ(IAXIS:KAXIS)=[(MESHES(NM).XF-MESHES(NM).XS)/MESHES(NM).IBAR ...
%                                  (MESHES(NM).YF-MESHES(NM).YS)/MESHES(NM).JBAR ...
%                                  (MESHES(NM).ZF-MESHES(NM).ZS)/MESHES(NM).KBAR];       
% end

return

