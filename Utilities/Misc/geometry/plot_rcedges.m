function [ierr]=plot_rcedges(NMESHES,MESHES,newfig)
IAXIS = 1; JAXIS = 2; KAXIS = 3;
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
    for IEDGE=1:MESHES(NM).IBM_NRCEDGE
        XCEN=MESHES(NM).IBM_RCEDGE(IEDGE).XCEN;
        DXI =MESHES(NM).IBM_RCEDGE(IEDGE).DXI;
        DX=0; DY=0; DZ=0;
        switch(MESHES(NM).IBM_RCEDGE(IEDGE).IJK(KAXIS+1))
            case(IAXIS)
                DX=DXI;
            case(JAXIS)
                DY=DXI; 
            case(KAXIS)
                DZ=DXI;
        end
        x=[XCEN(IAXIS)-DX/2 XCEN(IAXIS)+DX/2];
        y=[XCEN(JAXIS)-DY/2 XCEN(JAXIS)+DY/2];
        z=[XCEN(KAXIS)-DZ/2 XCEN(KAXIS)+DZ/2];

        plot3(x(1:2),y(1:2),z(1:2),colors(ct),'LineWidth',1)        
    end
end

if newfig
    axis equal; axis image;
    xlabel('X'); ylabel('Y'); zlabel('Z')
    view([45 45])
end
ierr=0;

return