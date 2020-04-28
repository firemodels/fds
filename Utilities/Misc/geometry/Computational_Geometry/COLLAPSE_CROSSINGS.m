function [ierr]=COLLAPSE_CROSSINGS(BODINT_PLANE2,X1AXIS,X2AXIS,X3AXIS,X3RAY,X1PLN,ititle)

global IAXIS JAXIS KAXIS GEOMEPS LOW_IND HIGH_IND
global IBM_GASPHASE IBM_SOLID IBM_UNDEFINED
global IBM_N_CRS IBM_MAXCROSS_X2
global IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_IS_CRS2_AUX IBM_BDNUM_CRS_AUX
global XAXIS


% Once all intersections and corresponding tags have been found, there
% might be points that lay on the same x2 location. Intersections type
% GG are dropped when other types are present at the same s. The remaining
% are reordered such that media continuity is preserved as the ray is
% covered for increasing s, by looking at the high side type of the
% adjacent intersection point to the left (if the intersection is the one
% with lowest s, the media to the left is type IBM_GASPHASE). Points of same
% type are collapsed. The final unique intersection type is obtained by
% using the LOW type of the first intersection and the HIGH type of the
% last intersection found at a given s.

IBM_N_CRS_AUX    = 0;
IBM_SVAR_CRS_AUX = 1./GEOMEPS*ones(1,IBM_MAXCROSS_X2); % svar = x2_intersection
IBM_IS_CRS2_AUX  = IBM_UNDEFINED*ones(HIGH_IND,IBM_MAXCROSS_X2); % Is the intersection an actual GS.
IBM_SEG_CRS_AUX  = zeros(1,IBM_MAXCROSS_X2);           % Segment containing the crossing.
IBM_SEG_TAN_AUX  = zeros(JAXIS,IBM_MAXCROSS_X2); % Segment orientation for each intersection.

% Count how many crossings with different SVAR:
CRS_NUM      = zeros(1,IBM_MAXCROSS_X2);
ICRS         = 1;
CRS_NUM(ICRS)= 1;
IND_CRS      = zeros(HIGH_IND,IBM_MAXCROSS_X2);
IND_CRS(LOW_IND, CRS_NUM(ICRS)) = ICRS-1;
IND_CRS(HIGH_IND,CRS_NUM(ICRS)) = IND_CRS(HIGH_IND,ICRS)+1;

for ICRS=2:IBM_N_CRS
   if ( abs(IBM_SVAR_CRS(ICRS)-IBM_SVAR_CRS(ICRS-1)) < GEOMEPS ) 
      CRS_NUM(ICRS) = CRS_NUM(ICRS-1);
   else
      CRS_NUM(ICRS) = CRS_NUM(ICRS-1)+1;
      IND_CRS(LOW_IND,CRS_NUM(ICRS)) = ICRS-1;
   end
   IND_CRS(HIGH_IND,CRS_NUM(ICRS)) = IND_CRS(HIGH_IND,CRS_NUM(ICRS))+1;
end

% Computation of IBM_BDNUM_CRS_AUX requires knowledge of how many different
% bodies reach an intersection:
BODNUM = zeros(1,IBM_MAXCROSS_X2);
for IDCR=1:CRS_NUM(IBM_N_CRS)
    
    % Load body numbers:
    for IDCR2=IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
       ISEG=IBM_SEG_CRS(IDCR2);
       if(ISEG > 0); BODNUM(IDCR2)=BODINT_PLANE2.INDSEG(4,ISEG); end
    end
    
    % Unique bodies:
    NUBD = 0;
    for IDCR2=IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
        if ( BODNUM(IDCR2)<1 ); continue; end
        if ((NUBD > 0) && any(UBOD(1:NUBD)==BODNUM(IDCR2))); continue; end        
        NUBD = NUBD + 1;
        UBOD(NUBD) = BODNUM(IDCR2);
    end
    % Now assign IBM_BDNUM_CRS_AUX(IDCR):
    SBOD = 0;
    for IUBD=1:NUBD
        % Drop extra intersections (same intersection type, same body):
        USE_INT_POINT(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)) = true;
        for ICRS1=IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
            if (~USE_INT_POINT(ICRS1)); continue; end % Don't use collapsed point as pivot.
            % Collapse GS or SG points:
            for ICRS2 = IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
                if ( (ICRS2==ICRS1) || ~USE_INT_POINT(ICRS2) ); continue; end % Don't use pivot, or collapsed point.
                if ( (IBM_IS_CRS2(LOW_IND ,ICRS1) == IBM_IS_CRS2(LOW_IND ,ICRS2)) && ...
                     (IBM_IS_CRS2(HIGH_IND,ICRS1) == IBM_IS_CRS2(HIGH_IND,ICRS2)) && ...
                      BODNUM(ICRS1) == BODNUM(ICRS2))
                    USE_INT_POINT(ICRS2) = false;
                end
            end
        end
        IBDNUM=0;
        for IDCR2=IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
            if(BODNUM(IDCR2) ~= UBOD(IUBD)); continue; end
            if(~USE_INT_POINT(IDCR2)); continue; end
            IBDNUM = IBDNUM + IBM_BDNUM_CRS(IDCR2);
        end
        if(IBDNUM ~= 0); SBOD = SBOD + sign(IBDNUM); end
    end

    if(IDCR == 1)
        IBM_BDNUM_CRS_AUX(IDCR) = SBOD;     
    else
        IBM_BDNUM_CRS_AUX(IDCR) = IBM_BDNUM_CRS_AUX(IDCR-1) + SBOD;
    end
end

% IDCR=1;
% IBM_BDNUM_CRS_AUX(IDCR) = max(0,max(IBM_BDNUM_CRS(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)))) ...
%                          +min(0,min(IBM_BDNUM_CRS(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR))));
% for IDCR=2:CRS_NUM(IBM_N_CRS)
%     IBM_BDNUM_CRS_AUX(IDCR) = IBM_BDNUM_CRS_AUX(IDCR-1) ...
%                              +max(0,max(IBM_BDNUM_CRS(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)))) ...
%                              +min(0,min(IBM_BDNUM_CRS(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR))));
% end


% This is where we merge intersections at same svar location (i.e. same CRS_NUM value):
% Loop over different crossings:
LEFT_MEDIA = IBM_GASPHASE; % Here we could change the initial LEFT_MEDIA to IBM_SOLID if needed. Would require adding
                           % IBM_BDNUM_CRS(LOW_IND,0) = 1, i.e crossed into SOLID at x2 -> -Inf.
for IDCR=1:CRS_NUM(IBM_N_CRS)

   IBM_N_CRS_AUX = IBM_N_CRS_AUX + 1;

   % Case of single crossing with new svar:
   if ( IND_CRS(HIGH_IND,IDCR) == 1 ) 

         ICRS =IND_CRS(LOW_IND,IDCR) + 1;

      if (ICRS > 1 && IBM_BDNUM_CRS_AUX(IDCR-1) > 0 && IBM_BDNUM_CRS_AUX(IDCR) > 0)  % Here test if we are already inside an Object.
         IBM_IS_CRS2(LOW_IND:HIGH_IND,ICRS) = IBM_SOLID;
      elseif ( IBM_IS_CRS2(LOW_IND,ICRS) ~= LEFT_MEDIA ) 

         % Check if this is a single point SGLS which was initially tagged as IBM_GASPHASE,
         % if so switch media type to LEFT_MEDIA
         if (IBM_SEG_CRS(ICRS) < 0) 
            IBM_IS_CRS2(LOW_IND:HIGH_IND,ICRS) = LEFT_MEDIA;
         else
            if(ititle==1)
            disp(['Error GET_X2INTERSECTIONS: IS_CRS(LOW_IND,ICRS) ~= LEFT_MEDIA, media continuity problem.'])
            disp(['XAXIS=' num2str(XAXIS) ', IBM_N_CRS=' num2str(IBM_N_CRS)])
            elseif(ititle==2)
            disp(['Error GET_IS_SOLID_PT: IS_CRS(LOW_IND,ICRS) ~= LEFT_MEDIA, media continuity problem.'])            
            disp(['XAXIS=' num2str(XAXIS) ', IBM_N_CRS=' num2str(IBM_N_CRS)])
            end
            disp(['X1AXIS,X1PLN=' num2str(X1AXIS) ', ' num2str(X1PLN) ...
                ', X2AXIS,X3AXIS=' num2str(X2AXIS) ', ' num2str(X3AXIS) ...
                ', RAY X3 POSITION=' num2str(X3RAY)])
            disp('IBM_SVAR_CRS, IBM_IS_CRS2')
            IBM_SVAR_CRS(1:IBM_N_CRS)
            IBM_IS_CRS2(1:2,1:IBM_N_CRS)
            disp('IBM_SVAR_CRS_AUX, IBM_IS_CRS2_AUX')
            IBM_SVAR_CRS_AUX(1:IBM_N_CRS_AUX)
            IBM_IS_CRS2_AUX(1:2,1:IBM_N_CRS_AUX)
            IBM_BDNUM_CRS_AUX(1:IBM_N_CRS_AUX)
            if(X2AXIS==IAXIS)
                x(1,:) = [-2 X3RAY X1PLN];
                x(2,:) = [ 2 X3RAY X1PLN];
            elseif(X2AXIS==JAXIS)
                x(1,:) = [X1PLN -2 X3RAY ];
                x(2,:) = [X1PLN  2 X3RAY ];
            else
                x(1,:) = [X3RAY X1PLN -2 ];
                x(2,:) = [X3RAY X1PLN  2 ];
            end
            plot3(x(:,IAXIS),x(:,JAXIS),x(:,KAXIS),'--k','LineWidth',2)
            pause
         end
      end

      IBM_SVAR_CRS_AUX(IBM_N_CRS_AUX)             = IBM_SVAR_CRS(ICRS);
      IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX)     = IBM_IS_CRS2(LOW_IND:HIGH_IND,ICRS);
      IBM_SEG_CRS_AUX(IBM_N_CRS_AUX)              = IBM_SEG_CRS(ICRS);
      IBM_SEG_TAN_AUX(IAXIS:JAXIS,IBM_N_CRS_AUX)  = IBM_SEG_TAN(IAXIS:JAXIS,ICRS);
      LEFT_MEDIA = IBM_IS_CRS2(HIGH_IND,ICRS);

      continue

   end
   
   % Case of several crossings with new svar:   
   % First test if any of the crossings lead to media type change, i.e
   % from IBM_SOLID to IBM_GASPHASE and viceversa.
   DROP_SS_GG = false;
   for ICRS=IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
      if ( IBM_IS_CRS2(LOW_IND,ICRS) ~= IBM_IS_CRS2(HIGH_IND,ICRS) ) 
         DROP_SS_GG = true;
         break
      end
   end

   % Variables related to new svar crossing:
   ICRS = IND_CRS(LOW_IND,IDCR) + 1;
   IBM_SVAR_CRS_AUX(IBM_N_CRS_AUX)             = IBM_SVAR_CRS(ICRS);
   IBM_SEG_CRS_AUX(IBM_N_CRS_AUX)              = IBM_SEG_CRS(ICRS);
   IBM_SEG_TAN_AUX(IAXIS:JAXIS,IBM_N_CRS_AUX)  = IBM_SEG_TAN(IAXIS:JAXIS,ICRS);
   
   % Case of intersection inside segment aligned with SVAR location, i.e.
   % intersection among two bodies or self intersection:
   ALGN_CROSS=false;
   for ICRS=IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
      if ( IBM_IS_CRS2(HIGH_IND+1,ICRS) ~= IBM_SOLID ); continue; end 
      IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX) = IBM_SOLID;
      ALGN_CROSS = true;
      break
   end
   if ( ALGN_CROSS ); continue; end

   % Now figure out the type of crossing:
   NOT_COUNTED = true;
   NCRS_REMAIN = IND_CRS(HIGH_IND,IDCR);
   if (DROP_SS_GG) 
 
      % Points of the same type are collapsed:
      USE_INT_POINT(IND_CRS(LOW_IND,IDCR)+1:IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)) = true;
      for ICRS1 = IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR) % Pivot Loop
         
         if (~USE_INT_POINT(ICRS1)); continue; end % Don't use collapsed point as pivot.

         % Collapse GS or SG points discriminating by body:
         for ICRS2 = IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
            if ( (ICRS2==ICRS1) || ~USE_INT_POINT(ICRS2) ); continue; end % Don't use pivot, or collapsed point.
            if ( (IBM_IS_CRS2(LOW_IND ,ICRS1) == IBM_IS_CRS2(LOW_IND ,ICRS2)) && ...
                 (IBM_IS_CRS2(HIGH_IND,ICRS1) == IBM_IS_CRS2(HIGH_IND,ICRS2)) && ...
                 (BODNUM(ICRS1) == BODNUM(ICRS2)) ) 
                USE_INT_POINT(ICRS2) = false;
            end
         end
         
      end
      
      % Left Side:
      FOUND_LEFT = false;
      IND_LEFT   = 0;
      IND_RIGHT  = 0;
      for ICRS=IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)
         % Case crossing type GG or SS, drop:
         if (IBM_IS_CRS2(LOW_IND,ICRS) == IBM_IS_CRS2(HIGH_IND,ICRS)); continue; end
         % Case collapsed point, drop:
         if (~USE_INT_POINT(ICRS)); continue; end

         IND_LEFT  =  IND_LEFT + IBM_IS_CRS2(LOW_IND,ICRS);
         IND_RIGHT = IND_RIGHT + IBM_IS_CRS2(HIGH_IND,ICRS);
      end
            
      
      if (IND_LEFT  ~= 0); IND_LEFT = sign(IND_LEFT); end
      if (IND_RIGHT ~= 0); IND_RIGHT = sign(IND_RIGHT); end
            
      if (IDCR>1 && IBM_BDNUM_CRS_AUX(IDCR-1) > 0 && IBM_BDNUM_CRS_AUX(IDCR) > 0)  % Test if we are inside an Object.
         IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX) = IBM_SOLID; % GS or SG.
      elseif (abs(IND_LEFT)+abs(IND_RIGHT) == 0)      % Same number of SG and GS crossings,
                                                      % both sides of the crossing
                                                      % defined as left_media:
         IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX)     = LEFT_MEDIA;
      elseif (IND_LEFT == LEFT_MEDIA) 
         IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX) = [ IND_LEFT, IND_RIGHT ]; % GS or SG.
      else
         if(ititle==1)
         disp('Error GET_X2INTERSECTIONS: DROP_SS_GG = true, Didn''t find left side continuity.')
         elseif(ititle==2)
         disp('Error GET_IS_SOLID_PT: DROP_SS_GG = true, Didn''t find left side continuity.')   
         disp(['XAXIS=' num2str(XAXIS)])
         end
         disp(['BODINT_PLANE2, NSGLS, NSEGS=',num2str(BODINT_PLANE2.NSGLS)  ', ' num2str(BODINT_PLANE2.NSEGS)])
         disp(['X1AXIS,X1PLN=' num2str(X1AXIS) ', ' num2str(X1PLN) ...
                ', X2AXIS,X3AXIS=' num2str(X2AXIS) ', ' num2str(X3AXIS) ...
                ', RAY X3 POSITION=' num2str(X3RAY)])
         disp('IBM_SVAR_CRS, IBM_IS_CRS2')
         IBM_SVAR_CRS(1:IBM_N_CRS)
         IBM_IS_CRS2(1:2,1:IBM_N_CRS)
         disp('IBM_SVAR_CRS_AUX, IBM_IS_CRS2_AUX')
         IBM_SVAR_CRS_AUX(1:IBM_N_CRS_AUX)
         IBM_IS_CRS2_AUX(1:2,1:IBM_N_CRS_AUX)
         pause
      end
      LEFT_MEDIA = IBM_IS_CRS2_AUX(HIGH_IND,IBM_N_CRS_AUX);

   else % Intersections are either GG or SS

      % Left side:
      FOUND_LEFT = false;
      for ICRS=IND_CRS(LOW_IND,IDCR)+1 : IND_CRS(LOW_IND,IDCR)+IND_CRS(HIGH_IND,IDCR)

         % Case GG or SS with IBM_IS_CRS2(LOW_IND,ICRS) == LEFT_MEDIA:
         % This collapses all types SS or GG that have the left side
         % type. Note they should all be one type (either GG or SS):
         if (IBM_IS_CRS2(LOW_IND,ICRS) == LEFT_MEDIA) 
            IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX) = IBM_IS_CRS2(LOW_IND:HIGH_IND,ICRS);
            NOT_COUNTED(ICRS) = false;
            NCRS_REMAIN = NCRS_REMAIN-1;
            FOUND_LEFT = true;
         end
      end

      if (IDCR>1 && IBM_BDNUM_CRS_AUX(IDCR-1) > 0 && IBM_BDNUM_CRS_AUX(IDCR) > 0)
         % Test if we are inside an Object.
         IBM_IS_CRS2_AUX(LOW_IND:HIGH_IND,IBM_N_CRS_AUX) = IBM_SOLID;
         LEFT_MEDIA = IBM_IS_CRS2_AUX(HIGH_IND,IBM_N_CRS_AUX);
         continue
      end
      
      if (~FOUND_LEFT)
          disp(['Error GET_X2INTERSECTIONS: DROP_SS_GG = false, Did not find left side continuity.'])
      end
      if (NCRS_REMAIN ~= 0) 
          disp(['Error GET_X2INTERSECTIONS: DROP_SS_GG = false, NCRS_REMAIN ~= 0.'])
      end
      LEFT_MEDIA = IBM_IS_CRS2_AUX(HIGH_IND,IBM_N_CRS_AUX);

   end

end

% Copy final results:
IBM_N_CRS    = IBM_N_CRS_AUX;
IBM_SVAR_CRS(1:IBM_MAXCROSS_X2)             = IBM_SVAR_CRS_AUX(1:IBM_MAXCROSS_X2);
IBM_SEG_CRS(1:IBM_MAXCROSS_X2)              = IBM_SEG_CRS_AUX(1:IBM_MAXCROSS_X2);
IBM_SEG_TAN(IAXIS:JAXIS,1:IBM_MAXCROSS_X2)  = IBM_SEG_TAN_AUX(IAXIS:JAXIS,1:IBM_MAXCROSS_X2);

for ICRS=1 : IBM_N_CRS
  IBM_IS_CRS(ICRS) = 2*( IBM_IS_CRS2_AUX(LOW_IND,ICRS) + 1 ) - IBM_IS_CRS2_AUX(HIGH_IND,ICRS);
end

ierr=0;

return