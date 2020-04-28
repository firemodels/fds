function [MESHES]=load_rcedges(basedir,casename,NMESHES,MESHES)
IAXIS = 1; JAXIS = 2; KAXIS = 3;

for NM=1:NMESHES
    DATA=load([basedir casename '_rcsegs_mesh_' num2str(NM,'%4.4d') '.dat']);
    if isempty(DATA)
        MESHES(NM).IBM_NRCEDGE=0;
        continue
    end
    MESHES(NM).IBM_NRCEDGE=length(DATA(:,1));
    for IEDGE=1:MESHES(NM).IBM_NRCEDGE
        MESHES(NM).IBM_RCEDGE(IEDGE).IJK(IAXIS:KAXIS+1)=DATA(IEDGE,IAXIS:KAXIS+1);
        MESHES(NM).IBM_RCEDGE(IEDGE).DXI               =DATA(IEDGE,5);
        MESHES(NM).IBM_RCEDGE(IEDGE).XCEN(IAXIS:KAXIS) =DATA(IEDGE,6:8);
    end
end

return