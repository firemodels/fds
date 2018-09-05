function [ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig)

IAXIS = 1; JAXIS = 2; KAXIS = 3;
ierr=1;

if newfig
    figure;
end
hold on
for IG=1:N_GEOMETRY
   % Geometry:
   [hg]=trisurf(GEOM(IG).WSELEM,GEOM(IG).XYZ(:,IAXIS),GEOM(IG).XYZ(:,JAXIS),GEOM(IG).XYZ(:,KAXIS));
end

if newfig
    axis equal; axis image;
    xlabel('X'); ylabel('Y'); zlabel('Z')
    view([45 45])
end
    
ierr=0;
return