function [MESHES,NMESHES]=load_meshes(basedir,casename)
IAXIS = 1; JAXIS = 2; KAXIS = 3;

DATA=load([basedir casename '_meshes.dat']);
NMESHES=length(DATA(:,1));
for NM=1:NMESHES
   MESHES(NM).IBAR = DATA(NM, 2); 
   MESHES(NM).JBAR = DATA(NM, 3); 
   MESHES(NM).KBAR = DATA(NM, 4);
   MESHES(NM).XS   = DATA(NM, 5);
   MESHES(NM).XF   = DATA(NM, 6);
   MESHES(NM).YS   = DATA(NM, 7);
   MESHES(NM).YF   = DATA(NM, 8);
   MESHES(NM).ZS   = DATA(NM, 9);
   MESHES(NM).ZF   = DATA(NM,10);
   MESHES(NM).DXYZ(IAXIS:KAXIS)=[(MESHES(NM).XF-MESHES(NM).XS)/MESHES(NM).IBAR ...
                                 (MESHES(NM).YF-MESHES(NM).YS)/MESHES(NM).JBAR ...
                                 (MESHES(NM).ZF-MESHES(NM).ZS)/MESHES(NM).KBAR];
       
end

return

