function [NFACE, CFELEM, BODTRI] = EAR_CLIP(NSEG,NVERT,SEG_CELL,XYZVERT,...
                                            INDIF,INDJF,INDKF,NPOLY,ILO_POLY,NSG_POLY)

global IAXIS JAXIS KAXIS
global NOD1 NOD2 NOD3
global GEOMEPS GEOM

pltflg = false; 

ICHK = 36;
JCHK = 17;
KCHK = 24;

if(pltflg && INDIF==ICHK && INDJF==JCHK && INDKF==KCHK)
    disp(['NPOLY=' num2str(NPOLY)])
    figure
    subplot(1,2,1)
    hold on
    a=0.0000001;
    colr = ['-or'; '-ob'; '-om';'-og';'-oc'];
    for ISEG=1:NSEG
        XYZSEG(:,NOD1)=XYZVERT(:,SEG_CELL(NOD1,ISEG));
        XYZSEG(:,NOD2)=XYZVERT(:,SEG_CELL(NOD2,ISEG));
        IG         =max(1,SEG_CELL(6,ISEG));
        plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),colr(IG,1:3),'LineWidth',2)
        text(0.5*(sum(XYZSEG(IAXIS,:)))+a,0.5*(sum(XYZSEG(JAXIS,:)))+a,...
             0.5*(sum(XYZSEG(KAXIS,:)))+a,num2str(SEG_CELL(4,ISEG)),'FontSize',10)
    end
    for IVERT=1:NVERT
        text(XYZVERT(IAXIS,IVERT)+a,XYZVERT(JAXIS,IVERT)+a,...
            XYZVERT(KAXIS,IVERT)+a,num2str(IVERT),'FontSize',10)
        
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %axis equal
    box on
    view([45 45])
    pause
end

% Compute segments director unit vectors and normals:
for ISEG=1:NSEG
    DV  = XYZVERT(IAXIS:KAXIS,SEG_CELL(NOD2,ISEG)) - ...
          XYZVERT(IAXIS:KAXIS,SEG_CELL(NOD1,ISEG));
    LEN(ISEG) = norm(DV);
    N(IAXIS:KAXIS,ISEG) = 1/LEN(ISEG) * DV;
end

% First sweep across segments defining triangles for all consecutive
% segments with same triangle and body: 
SEG_CELL2 = zeros(6,NSEG);
NFACE     = 0;
CFELEM    = [];
BODTRI    = [];
% Ear clipping algorithm by polyline:
for IPOLY=1:NPOLY
    
    ILO       = ILO_POLY(IPOLY)+1;
    NSGP      = NSG_POLY(IPOLY);
    IHI       = ILO_POLY(IPOLY)+NSGP;
    SEG_CELL2(:,1:NSGP) = SEG_CELL(:,ILO:IHI);
    SEG_FLAG  = zeros(1,NSGP);
    LEFT_SEGS = NSGP; 
    
    for COUNTEXT=1:3  % Search segmets first that belong to same triangle (1),
        %               second that belong to same body (2),
        %                third all the rest.
        for COUNT=1:2 % Search first last uncounted segment (1), second the rest.
            if(LEFT_SEGS < 3); break; end % should break out of COUNTEXT loop.
            if (COUNT==1)
                ISEG  = NSGP-1;
                for ISEG1=1:NSGP
                    if (~SEG_FLAG(ISEG1)); break; end
                end
            else
                ISEG  = 0;
            end
            while (ISEG < NSGP)
                ISEG = ISEG + 1;
                
                if (SEG_FLAG(ISEG)); continue; end
                
                FOUND_ISEG1 = 0;
                if (COUNT==1)
                    if (~SEG_FLAG(ISEG1)); FOUND_ISEG1 = 1; end
                else
                    for ISEG1=ISEG+1:NSGP
                        if (~SEG_FLAG(ISEG1))
                            FOUND_ISEG1 = 1;
                            break
                        end
                    end
                end
                if(~FOUND_ISEG1); continue; end
                
                TRI   = 0;
                % Test if triangle given by ISEG ISEG+1 DIAG is valid.
                % First, drop if Body not the same:
                if(COUNTEXT<3 && (SEG_CELL2(6,ISEG) ~= SEG_CELL2(6,ISEG1))); continue; end
                
                % Second, drop if segments are on the same line:
                if(abs(abs(dot(N(IAXIS:KAXIS,ISEG),N(IAXIS:KAXIS,ISEG1)))-1.) < 1.e-12); continue; end
                
                % Now drop if triangles don't match:
                TWOTRI = 0;
                if(COUNTEXT<3)
                    if (SEG_CELL2(4,ISEG) ~=0 && (SEG_CELL2(4,ISEG) == SEG_CELL2(4,ISEG1) || ...
                            SEG_CELL2(4,ISEG) == SEG_CELL2(5,ISEG1)) )
                        TWOTRI=1;
                        TRI = SEG_CELL2(4,ISEG);
                        BOD = SEG_CELL2(6,ISEG);
                    elseif (SEG_CELL2(5,ISEG) ~=0 && (SEG_CELL2(5,ISEG) == SEG_CELL2(4,ISEG1) || ...
                            SEG_CELL2(5,ISEG) == SEG_CELL2(5,ISEG1)))
                        TWOTRI=1;
                        TRI = SEG_CELL2(5,ISEG);
                        BOD = SEG_CELL2(6,ISEG);
                    end
                end
                if(COUNTEXT~=1 && TRI==0)
                    % Define TRI as the longest seg one:
                    if ( LEN(ISEG) >= LEN(ISEG1) )
                        TRI = SEG_CELL2(4,ISEG);
                        BOD = SEG_CELL2(6,ISEG);
                    else
                        TRI = SEG_CELL2(4,ISEG1);
                        BOD = SEG_CELL2(6,ISEG1);
                    end
                end
                
                if(TRI == 0 )
                    continue
                else % Found two segments with matching triangle.
                    
                    % Test that triangle found is not internal to GEOMs:
                    CONN = [SEG_CELL2(1:2,ISEG)' SEG_CELL2(2,ISEG1)];              
                    if(TWOTRI)
                        NP = GEOM(BOD).FACES_NORMAL(IAXIS:KAXIS,TRI)';    
                        XP(IAXIS:KAXIS) = ...
                        1/3*( XYZVERT(IAXIS:KAXIS,CONN(NOD1)) + ...
                              XYZVERT(IAXIS:KAXIS,CONN(NOD2)) + ...
                              XYZVERT(IAXIS:KAXIS,CONN(NOD3)))'+ 10.*GEOMEPS*NP;
                        [val,XAXIS] = max(abs(NP(IAXIS:KAXIS)));
                        [IS_SOLID]=GET_IS_SOLID_3D(XAXIS,XP,INDIF,INDJF,INDKF);
                        if(IS_SOLID); continue; end
                    end
                    
                    NFACE = NFACE + 1;
                    CFELEM(1:4,NFACE) = [3 CONN]';
                    BODTRI(:,NFACE)   = [BOD TRI]';
                    SEG_CELL2(:,ISEG) = [SEG_CELL2(1,ISEG) SEG_CELL2(2,ISEG1) 1 TRI 0 BOD]';
                    
                    XYZ1=XYZVERT(IAXIS:KAXIS,SEG_CELL2(1,ISEG));
                    XYZ2=XYZVERT(IAXIS:KAXIS,SEG_CELL2(2,ISEG));
                    DV  = XYZ2 - XYZ1;
                    LEN(ISEG) = norm(DV);
                    N(IAXIS:KAXIS,ISEG) = 1/LEN(ISEG) * DV;
                    
                    % Erase Segment from
                    SEG_CELL2(:,ISEG1)  = 0;
                    SEG_FLAG(ISEG1)     = 1;
                    N(IAXIS:KAXIS,ISEG1)= 0.;
                    LEFT_SEGS = LEFT_SEGS - 1;
                    
                    if(COUNT~=1); ISEG = ISEG - 1; end
                end
            end
        end
    end

end

if(pltflg && INDIF==ICHK && INDJF==JCHK && INDKF==KCHK)
    subplot(1,2,2)
    hold on
    for ISEG=1:NSGP
        if(SEG_FLAG(ISEG)); continue; end
        XYZSEG(:,NOD1)=XYZVERT(:,SEG_CELL2(NOD1,ISEG));
        XYZSEG(:,NOD2)=XYZVERT(:,SEG_CELL2(NOD2,ISEG));
        IG         =max(1,SEG_CELL2(6,ISEG));
        plot3(XYZSEG(IAXIS,:),XYZSEG(JAXIS,:),XYZSEG(KAXIS,:),'--k','LineWidth',2)
        text(0.5*(sum(XYZSEG(IAXIS,:)))+a,0.5*(sum(XYZSEG(JAXIS,:)))+a,...
             0.5*(sum(XYZSEG(KAXIS,:)))+a,num2str(SEG_CELL(4,ISEG)),'FontSize',10)
    end
    % Plot faces:
    for ICF=1:NFACE
        XYZFC(:,NOD1)=XYZVERT(:,CFELEM(2,ICF));
        XYZFC(:,NOD2)=XYZVERT(:,CFELEM(3,ICF));
        XYZFC(:,NOD3)=XYZVERT(:,CFELEM(4,ICF));
        IG         =BODTRI(1,ICF);
        [hp]=patch(XYZFC(IAXIS,:),XYZFC(JAXIS,:),XYZFC(KAXIS,:),colr(IG,1:3));
        set(hp,'FaceAlpha',0.1)
        text(1/3*(sum(XYZFC(IAXIS,:)))+a,1/3*(sum(XYZFC(JAXIS,:)))+a,...
             1/3*(sum(XYZFC(KAXIS,:)))+a,num2str(BODTRI(2,ICF)),'FontSize',10)
        pause
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %axis equal
    box on
    view([45 45])
    pause
end

return