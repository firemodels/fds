function [ierr]=plot_geometries(N_GEOMETRY,GEOM,newfig)

IAXIS = 1; JAXIS = 2; KAXIS = 3;
NOD1=1; NOD2=2; NOD3=3;
ierr=1;

if newfig
    figure;
end
hold on
for IG=1:N_GEOMETRY
   % Geometry:
   [hg]=trisurf(GEOM(IG).WSELEM,GEOM(IG).XYZ(:,IAXIS),GEOM(IG).XYZ(:,JAXIS),GEOM(IG).XYZ(:,KAXIS));
   set(hg,'FaceAlpha',0.3)
   V=GEOM(IG).XYZ; F=GEOM(IG).WSELEM; N=GEOM(IG).FACES_NORMAL';
   nfaces=GEOM(IG).N_FACES;
   a=0.05;
   xyz1=1/3*(V(F(1:nfaces,NOD1),IAXIS:KAXIS) + ...
             V(F(1:nfaces,NOD2),IAXIS:KAXIS) + ...
             V(F(1:nfaces,NOD3),IAXIS:KAXIS));
   xyz2=xyz1(1:nfaces,IAXIS:KAXIS)+a*N(1:nfaces,IAXIS:KAXIS);     
   plot3([xyz1(1:nfaces,IAXIS) xyz2(1:nfaces,IAXIS)]', ...
         [xyz1(1:nfaces,JAXIS) xyz2(1:nfaces,JAXIS)]', ...
         [xyz1(1:nfaces,KAXIS) xyz2(1:nfaces,KAXIS)]','r'); 
   
end

if newfig
    axis equal; axis image;
    xlabel('X'); ylabel('Y'); zlabel('Z')
    view([45 45])
end
    
ierr=0;
return