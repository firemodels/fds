function [GEOM,N_GEOMETRY]=load_geometries(basedir,casename);

N_GEOMETRY=load([basedir casename '_num_geometries.dat']);
for IG=1:N_GEOMETRY
   GEOM(IG).XYZ   =load([basedir casename '_geometry_' num2str(IG,'%4.4d') '_verts.dat']); 
   GEOM(IG).WSELEM=load([basedir casename '_geometry_' num2str(IG,'%4.4d') '_faces.dat']); 
end

return