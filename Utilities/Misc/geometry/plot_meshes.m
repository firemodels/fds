function [ierr]=plot_meshes(NMESHES,MESHES,newfig)

ierr=1;
colors=['b'; 'g'; 'r'; 'c'; 'm'; 'y'; 'k'];
NT=ceil(NMESHES/7);
for i=1:NT-1
    colors = [colors; colors];
end

if newfig
    figure;
end
hold on
ct=0;
for NM=1:NMESHES
   
   ct=ct+1;

   % Create x, y, z lists for mesh edges:
   x( 1,1:2)=[MESHES(NM).XS MESHES(NM).XF];
   y( 1,1:2)=[MESHES(NM).YS MESHES(NM).YS];
   z( 1,1:2)=[MESHES(NM).ZS MESHES(NM).ZS];
   
   x( 2,1:2)=[MESHES(NM).XS MESHES(NM).XF];
   y( 2,1:2)=[MESHES(NM).YF MESHES(NM).YF];
   z( 2,1:2)=[MESHES(NM).ZS MESHES(NM).ZS];

   x( 3,1:2)=[MESHES(NM).XS MESHES(NM).XF];
   y( 3,1:2)=[MESHES(NM).YS MESHES(NM).YS];
   z( 3,1:2)=[MESHES(NM).ZF MESHES(NM).ZF];
   
   x( 4,1:2)=[MESHES(NM).XS MESHES(NM).XF];
   y( 4,1:2)=[MESHES(NM).YF MESHES(NM).YF];
   z( 4,1:2)=[MESHES(NM).ZF MESHES(NM).ZF];
   
   x( 5,1:2)=[MESHES(NM).XS MESHES(NM).XS];
   y( 5,1:2)=[MESHES(NM).YS MESHES(NM).YF];
   z( 5,1:2)=[MESHES(NM).ZS MESHES(NM).ZS];
   
   x( 6,1:2)=[MESHES(NM).XF MESHES(NM).XF];
   y( 6,1:2)=[MESHES(NM).YS MESHES(NM).YF];
   z( 6,1:2)=[MESHES(NM).ZS MESHES(NM).ZS];
   
   x( 7,1:2)=[MESHES(NM).XS MESHES(NM).XS];
   y( 7,1:2)=[MESHES(NM).YS MESHES(NM).YF];
   z( 7,1:2)=[MESHES(NM).ZF MESHES(NM).ZF];
   
   x( 8,1:2)=[MESHES(NM).XF MESHES(NM).XF];
   y( 8,1:2)=[MESHES(NM).YS MESHES(NM).YF];
   z( 8,1:2)=[MESHES(NM).ZF MESHES(NM).ZF];
   
   x( 9,1:2)=[MESHES(NM).XS MESHES(NM).XS];
   y( 9,1:2)=[MESHES(NM).YS MESHES(NM).YS];
   z( 9,1:2)=[MESHES(NM).ZS MESHES(NM).ZF];
   
   x(10,1:2)=[MESHES(NM).XF MESHES(NM).XF];
   y(10,1:2)=[MESHES(NM).YS MESHES(NM).YS];
   z(10,1:2)=[MESHES(NM).ZS MESHES(NM).ZF];
   
   x(11,1:2)=[MESHES(NM).XS MESHES(NM).XS];
   y(11,1:2)=[MESHES(NM).YF MESHES(NM).YF];
   z(11,1:2)=[MESHES(NM).ZS MESHES(NM).ZF];
   
   x(12,1:2)=[MESHES(NM).XF MESHES(NM).XF];
   y(12,1:2)=[MESHES(NM).YF MESHES(NM).YF];
   z(12,1:2)=[MESHES(NM).ZS MESHES(NM).ZF];
   
   for iseg=1:12
      plot3(x(iseg,1:2),y(iseg,1:2),z(iseg,1:2),colors(ct),'LineWidth',2)
   end

end

if newfig
    axis equal; axis image;
    xlabel('X'); ylabel('Y'); zlabel('Z')
    view([45 45])
end
    

ierr=0;
return