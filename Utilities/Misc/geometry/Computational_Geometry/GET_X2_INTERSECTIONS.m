function [ierr] = GET_X2_INTERSECTIONS(X1AXIS,X2AXIS,X3AXIS,X3RAY,X1PLN)

global IAXIS JAXIS KAXIS NOD1 NOD2 GEOMEPS LOW_IND HIGH_IND
global IBM_GASPHASE IBM_SOLID IBM_UNDEFINED
global IBM_N_CRS IBM_MAXCROSS_X2
global BODINT_PLANE
global IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_IS_CRS2_AUX IBM_BDNUM_CRS_AUX
global IBM_OR_CRS
ierr = 1;

% Initialize crossings arrays:
IBM_N_CRS    = 0;
IBM_SVAR_CRS = 1./GEOMEPS*ones(1,IBM_MAXCROSS_X2);
IBM_IS_CRS   = IBM_UNDEFINED*ones(1,IBM_MAXCROSS_X2);
IBM_IS_CRS2  = IBM_UNDEFINED*ones(3,IBM_MAXCROSS_X2);
IBM_SEG_TAN  = zeros(JAXIS,IBM_MAXCROSS_X2);
IBM_SEG_CRS  = zeros(1,IBM_MAXCROSS_X2);
IBM_BDNUM_CRS= zeros(LOW_IND,IBM_MAXCROSS_X2);
IBM_BDNUM_CRS_AUX= zeros(LOW_IND,IBM_MAXCROSS_X2);

% First Single points:
% Treat them as [GASPHASE GASPHASE] crossings:
for ISGL=1:BODINT_PLANE.NSGLS

   SGL = BODINT_PLANE.SGLS(NOD1,ISGL);
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SGL);

   % x2-x3 coordinates of point:
   X2_1 = XYZ1(X2AXIS);
   X3_1 = XYZ1(X3AXIS);

   % Dot product dot(X_1-XRAY,e3)
   DOT1 = X3_1-X3RAY;
   if (abs(DOT1) <= GEOMEPS); DOT1=0.; end
   if ( abs(DOT1) == 0. )
       % Point 1:
       SVARI = X2_1;
       ICRSI(LOW_IND:HIGH_IND+1) = [IBM_GASPHASE, IBM_GASPHASE, IBM_UNDEFINED]; % Both IBM_GASPHASE for SGL point.
       SCRSI = -ISGL;
       STANI(IAXIS:JAXIS)  = 0.;

       % Insertion sort:
       [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI); % Modifies crossings arrays.
   end   
end

if(BODINT_PLANE.NSEGS > 0)
DELBIN = BODINT_PLANE.TBAXIS(X3AXIS).DELBIN;
AVAL   = (X3RAY-GEOMEPS-BODINT_PLANE.BOX(LOW_IND,X3AXIS))/DELBIN;
ILO_BIN= max(1,ceil(sign(AVAL)*min(2*BODINT_PLANE.TBAXIS(X3AXIS).N_BINS,abs(AVAL))));
AVAL   = (X3RAY+GEOMEPS-BODINT_PLANE.BOX(LOW_IND,X3AXIS))/DELBIN;
IHI_BIN= min(BODINT_PLANE.TBAXIS(X3AXIS).N_BINS, ...
         ceil(sign(AVAL)*min(2*BODINT_PLANE.TBAXIS(X3AXIS).N_BINS,abs(AVAL))));
for IBIN=ILO_BIN:IHI_BIN

   if (X3RAY < BODINT_PLANE.TBAXIS(X3AXIS).TRIBIN(IBIN).X_LOW-GEOMEPS); continue; end
   if (X3RAY > BODINT_PLANE.TBAXIS(X3AXIS).TRIBIN(IBIN).X_HIGH+GEOMEPS); continue; end

   for IISEG=1:BODINT_PLANE.TBAXIS(X3AXIS).TRIBIN(IBIN).NTL

   ISEG = BODINT_PLANE.TBAXIS(X3AXIS).TRIBIN(IBIN).TRI_LIST(IISEG);
%for ISEG=1:BODINT_PLANE.NSEGS
    
   SEG(NOD1:NOD2)    = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD1));
   XYZ2(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD2));

   % x2,x3 coordinates of segment:
   X2_1 = XYZ1(X2AXIS);
   X3_1 = XYZ1(X3AXIS); % Lower index endpoint.
   X2_2 = XYZ2(X2AXIS);
   X3_2 = XYZ2(X3AXIS); % Upper index endpoint.

   % Is segment aligned with x3 direction?
   %BODINT_PLANE.X3ALIGNED(ISEG) = (abs(X2_2-X2_1) < GEOMEPS);
   % Is segment aligned with x2 rays?:
   %BODINT_PLANE.X2ALIGNED(ISEG) = (abs(X3_2-X3_1) < GEOMEPS);

   % First Test if the whole segment is on one side of the Ray:
   % Test segment crosses the ray, or is in geomepsilon proximity
   % of it:
   X3MIN = min(X3_1,X3_2);
   X3MAX = max(X3_1,X3_2);
   OUTRAY=(((X3RAY-X3MAX) > GEOMEPS) || ((X3MIN-X3RAY) > GEOMEPS));

   if (OUTRAY); continue; end

   DOT1 = X3_1-X3RAY;
   DOT2 = X3_2-X3RAY;

   if (abs(DOT1) <= GEOMEPS); DOT1 = 0.; end
   if (abs(DOT2) <= GEOMEPS); DOT2 = 0.; end

   % Segment tangent unit vector.
   DV12(IAXIS:JAXIS) = XYZ2( [ X2AXIS, X3AXIS ] ) - XYZ1( [ X2AXIS, X3AXIS ] );
   MODTI = sqrt( DV12(IAXIS)^2. + DV12(JAXIS)^2. );
   STANI(IAXIS:JAXIS)  = DV12(IAXIS:JAXIS) * MODTI^(-1.);
   NOMLI(IAXIS:JAXIS)  = [ STANI(JAXIS), -STANI(IAXIS) ];
   ISSEG(LOW_IND:HIGH_IND) = BODINT_PLANE.SEGTYPE(LOW_IND:HIGH_IND,ISEG);

   % For x2, in local x2-x3 coords e2=(1,0):
   GAM(LOW_IND) = (1 + sign(NOMLI(IAXIS))) / 2; % (1+SIGN(DOT_PRODUCT(NOMLI,e2)))/2;
   GAM(HIGH_IND)= (1 - sign(NOMLI(IAXIS))) / 2; % (1-SIGN(DOT_PRODUCT(NOMLI,e2)))/2;

   % Test if whole segment is in ray, if so add segment nodes as crossings:
   if ( (abs(DOT1)+abs(DOT2)) == 0. )

      % Count both points as crossings:
      % Point 1:
      SVARI = min(X2_1,X2_2);
      ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_GASPHASE, IBM_SOLID, IBM_UNDEFINED];
      SCRSI = ISEG;

      % Insertion sort:
      [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);

      
      for ICR=2:BODINT_PLANE.NBCROSS(ISEG)-1
          SVARI = X2_1 + BODINT_PLANE.SVAR(ICR,ISEG)*STANI(IAXIS);
          ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_SOLID, IBM_SOLID, IBM_SOLID];
          SCRSI = ISEG;
          % Insertion sort:
          [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
      end
      
      
      % Point 2:
      SVARI = max(X2_1,X2_2);
      ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_SOLID, IBM_GASPHASE, IBM_UNDEFINED ];
      SCRSI = ISEG;

      % Insertion sort:
      [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);

      continue

   end

   % Now nodes individually:
   if ( abs(DOT1) == 0. ) 

      % Point 1:
      SVARI = X2_1;

      % LOW and HIGH media type, using the segment definition:
      ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
      ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
      ICRSI(HIGH_IND+1) = IBM_UNDEFINED;
      SCRSI = ISEG;

      % Insertion sort:
      [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);

      continue

   end
   if ( abs(DOT2) == 0. ) 

      % Point 2:
      SVARI = X2_2;

      % LOW and HIGH_IND media type, using the segment definition:
      ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
      ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
      ICRSI(HIGH_IND+1) = IBM_UNDEFINED;
      SCRSI = ISEG;

      % Insertion sort:
      [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);

      continue

   end

   % Finally regular case:
   % Points 1 on one side of ray, point 2 on the other:
   % if ((DOT1 > 0. .AND. DOT2 < 0.) || (DOT1 < 0. .AND. DOT2 > 0.))
   if ( DOT1*DOT2 < 0. ) 

      % Intersection Point along segment:
      SVARI = X2_1 + (X3RAY-X3_1) * (X2_2-X2_1) / (X3_2-X3_1);

      % LOW and HIGH media type, using the segment definition:
      ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
      ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
      ICRSI(HIGH_IND+1) = IBM_UNDEFINED;
      SCRSI = ISEG;

      % Insertion sort:
      [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);

      continue

   end

   disp(['Error GET_X2INTERSECTIONS: Missed segment=' num2str(ISEG)])

   end
end

end

% Do we have any intersections?
if ( IBM_N_CRS == 0 ); return; end

% Collapse crossings to single SVARs:
[ierr2]=COLLAPSE_CROSSINGS(BODINT_PLANE,X1AXIS,X2AXIS,X3AXIS,X3RAY,X1PLN,1);

ierr = 0;

return