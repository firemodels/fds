function [IS_SOLID] = GET_IS_SOLID_PT(BODINT_PLANE2,X1AXIS,X2AXIS,X3AXIS,XY,NVEC,X1PLN)


global IAXIS JAXIS KAXIS NOD1 NOD2 GEOMEPS LOW_IND HIGH_IND
global IBM_GASPHASE IBM_SOLID IBM_UNDEFINED
global IBM_N_CRS IBM_MAXCROSS_X2
global IBM_SVAR_CRS IBM_IS_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS IBM_IS_CRS2_AUX IBM_BDNUM_CRS_AUX
global XAXIS
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

% Define crossings
if(abs(NVEC(IAXIS)) > abs(NVEC(JAXIS))) % Do X2 ray
   SCEN=XY(IAXIS);
   XRAY=XY(JAXIS);
   XAXIS=X3AXIS;

   DELBIN = BODINT_PLANE2.TBAXIS(XAXIS).DELBIN;
   AVAL   = (XRAY-GEOMEPS-BODINT_PLANE2.BOX(LOW_IND,XAXIS))/DELBIN;
   ILO_BIN= max(1,ceil(sign(AVAL)*min(2*BODINT_PLANE2.TBAXIS(XAXIS).N_BINS,abs(AVAL))));
   AVAL   = (XRAY+GEOMEPS-BODINT_PLANE2.BOX(LOW_IND,XAXIS))/DELBIN;
   IHI_BIN= min(BODINT_PLANE2.TBAXIS(XAXIS).N_BINS, ...
            ceil(sign(AVAL)*min(2*BODINT_PLANE2.TBAXIS(XAXIS).N_BINS,abs(AVAL))));
   for IBIN=ILO_BIN:IHI_BIN
      if (XRAY < BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).X_LOW-GEOMEPS); continue; end
      if (XRAY > BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).X_HIGH+GEOMEPS); continue; end
      for IISEG=1:BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).NTL
         ISEG = BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).TRI_LIST(IISEG);
         SEG(NOD1:NOD2)    = BODINT_PLANE2.SEGS(NOD1:NOD2,ISEG);
         XYZ1(IAXIS:KAXIS) = BODINT_PLANE2.XYZ(IAXIS:KAXIS,SEG(NOD1));
         XYZ2(IAXIS:KAXIS) = BODINT_PLANE2.XYZ(IAXIS:KAXIS,SEG(NOD2));

         % x2,x3 coordinates of segment:
         X2_1 = XYZ1(X2AXIS);
         X3_1 = XYZ1(X3AXIS); % Lower index endpoint.
         X2_2 = XYZ2(X2AXIS);
         X3_2 = XYZ2(X3AXIS); % Upper index endpoint.

         % First Test if the whole segment is on one side of the Ray:
         % Test segment crosses the ray, or is in geomepsilon proximity
         % of it:
         X3MIN = min(X3_1,X3_2);
         X3MAX = max(X3_1,X3_2);
         OUTRAY=(((XRAY-X3MAX) > GEOMEPS) || ((X3MIN-XRAY) > GEOMEPS));

         if (OUTRAY); continue; end
         DOT1 = X3_1-XRAY;          DOT2 = X3_2-XRAY;
         if (abs(DOT1) <= GEOMEPS); DOT1 = 0.; end
         if (abs(DOT2) <= GEOMEPS); DOT2 = 0.; end

         % Segment tangent unit vector.
         DV12(IAXIS:JAXIS) = XYZ2( [ X2AXIS, X3AXIS ] ) - XYZ1( [ X2AXIS, X3AXIS ] );
         MODTI = sqrt( DV12(IAXIS)^2. + DV12(JAXIS)^2. );
         STANI(IAXIS:JAXIS)  = DV12(IAXIS:JAXIS) * MODTI^(-1.);
         NOMLI(IAXIS:JAXIS)  = [ STANI(JAXIS), -STANI(IAXIS) ];
         ISSEG(LOW_IND:HIGH_IND) = BODINT_PLANE2.SEGTYPE(LOW_IND:HIGH_IND,ISEG);

         % For x2, in local x2-x3 coords e2=(1,0):
         GAM(LOW_IND) = (1 + sign(NOMLI(IAXIS))) / 2; % (1+sign(DOT_PRODUCT(NOMLI,e2)))/2;
         GAM(HIGH_IND)= (1 - sign(NOMLI(IAXIS))) / 2; % (1-sign(DOT_PRODUCT(NOMLI,e2)))/2;

         % Test if whole segment is in ray, if so add segment nodes as crossings:
         if ( (abs(DOT1)+abs(DOT2)) == 0. )
            % Count both points as crossings:
            % Point 1:
            SVARI = min(X2_1,X2_2);
            ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_GASPHASE, IBM_SOLID, IBM_UNDEFINED ];
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            for ICR=2:BODINT_PLANE2.NBCROSS(ISEG)-1
                SVARI = X2_1 + BODINT_PLANE2.SVAR(ICR,ISEG)*STANI(IAXIS);
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
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
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
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         % Finally regular case:
         % Points 1 on one side of ray, point 2 on the other:
         if ( DOT1*DOT2 < 0. ) 
            % Intersection Point along segment:
            % DS    = (XRAY-X3_1) / (X3_2-X3_1)
            % SVARI = X2_1 + DS*(X2_2-X2_1)
            SVARI = X2_1 + (XRAY-X3_1) * (X2_2-X2_1) / (X3_2-X3_1);
            % LOW and HIGH media type, using the segment definition:
            ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
            ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         disp(['Error GET_IS_SOLID_PT NVEC(IAXIS): Missed segment=' num2str(ISEG)])
      end
    end 
    
else % Do X3 ray
   SCEN=XY(JAXIS);
   XRAY=XY(IAXIS);
   XAXIS=X2AXIS;
    
   DELBIN = BODINT_PLANE2.TBAXIS(XAXIS).DELBIN;
   AVAL   = (XRAY-GEOMEPS-BODINT_PLANE2.BOX(LOW_IND,XAXIS))/DELBIN;
   ILO_BIN= max(1,ceil(sign(AVAL)*min(2*BODINT_PLANE2.TBAXIS(XAXIS).N_BINS,abs(AVAL))));
   AVAL   = (XRAY+GEOMEPS-BODINT_PLANE2.BOX(LOW_IND,XAXIS))/DELBIN;
   IHI_BIN= min(BODINT_PLANE2.TBAXIS(XAXIS).N_BINS, ...
            ceil(sign(AVAL)*min(2*BODINT_PLANE2.TBAXIS(XAXIS).N_BINS,abs(AVAL))));
   for IBIN=ILO_BIN:IHI_BIN
      if (XRAY < BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).X_LOW-GEOMEPS); continue; end
      if (XRAY > BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).X_HIGH+GEOMEPS); continue; end
      for IISEG=1:BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).NTL
         ISEG = BODINT_PLANE2.TBAXIS(XAXIS).TRIBIN(IBIN).TRI_LIST(IISEG);
         SEG(NOD1:NOD2)    = BODINT_PLANE2.SEGS(NOD1:NOD2,ISEG);
         XYZ1(IAXIS:KAXIS) = BODINT_PLANE2.XYZ(IAXIS:KAXIS,SEG(NOD1));
         XYZ2(IAXIS:KAXIS) = BODINT_PLANE2.XYZ(IAXIS:KAXIS,SEG(NOD2));

         % x2,x3 coordinates of segment:
         X2_1 = XYZ1(X2AXIS);
         X3_1 = XYZ1(X3AXIS); % Lower index endpoint.
         X2_2 = XYZ2(X2AXIS);
         X3_2 = XYZ2(X3AXIS); % Upper index endpoint.

         % First Test if the whole segment is on one side of the Ray:
         % Test segment crosses the ray, or is in geomepsilon proximity
         % of it:
         X2MIN = min(X2_1,X2_2);
         X2MAX = max(X2_1,X2_2);
         OUTRAY=(((XRAY-X2MAX) > GEOMEPS) || ((X2MIN-XRAY) > GEOMEPS));

         if (OUTRAY); continue; end
         DOT1 = X2_1-XRAY;          DOT2 = X2_2-XRAY;
         if (abs(DOT1) <= GEOMEPS); DOT1 = 0.; end
         if (abs(DOT2) <= GEOMEPS); DOT2 = 0.; end

         % Segment tangent unit vector.
         DV12(IAXIS:JAXIS) = XYZ2( [ X2AXIS, X3AXIS ] ) - XYZ1( [ X2AXIS, X3AXIS ] );
         MODTI = sqrt( DV12(IAXIS)^2. + DV12(JAXIS)^2. );
         STANI(IAXIS:JAXIS)  = DV12(IAXIS:JAXIS) * MODTI^(-1.);
         NOMLI(IAXIS:JAXIS)  = [ STANI(JAXIS), -STANI(IAXIS) ];
         ISSEG(LOW_IND:HIGH_IND) = BODINT_PLANE2.SEGTYPE(LOW_IND:HIGH_IND,ISEG);         
         
         % For x3, in local x2-x3 coords e2=(0,1):
         GAM(LOW_IND) = (1 + sign(NOMLI(JAXIS))) / 2; % (1+sign(DOT_PRODUCT(NOMLI,e2)))/2;
         GAM(HIGH_IND)= (1 - sign(NOMLI(JAXIS))) / 2; % (1-sign(DOT_PRODUCT(NOMLI,e2)))/2;

         % Test if whole segment is in ray, if so add segment nodes as crossings:
         if ( (abs(DOT1)+abs(DOT2)) == 0. )
            % Count both points as crossings:
            % Point 1:
            SVARI = min(X3_1,X3_2);
            ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_GASPHASE, IBM_SOLID, IBM_UNDEFINED ];
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            for ICR=2:BODINT_PLANE2.NBCROSS(ISEG)-1
                SVARI = X3_1 + BODINT_PLANE2.SVAR(ICR,ISEG)*STANI(JAXIS);;
                ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_SOLID, IBM_SOLID, IBM_SOLID];
                SCRSI = ISEG;
                % Insertion sort:
                [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            end
            % Point 2:
            SVARI = max(X3_1,X3_2);
            ICRSI(LOW_IND:HIGH_IND+1) = [ IBM_SOLID, IBM_GASPHASE, IBM_UNDEFINED ];
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         % Now nodes individually:
         if ( abs(DOT1) == 0. ) 
            % Point 1:
            SVARI = X3_1;
            % LOW and HIGH media type, using the segment definition:
            ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
            ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         if ( abs(DOT2) == 0. ) 
            % Point 2:
            SVARI = X3_2;
            % LOW and HIGH_IND media type, using the segment definition:
            ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
            ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         % Finally regular case:
         % Points 1 on one side of ray, point 2 on the other:
         if ( DOT1*DOT2 < 0. ) 
            % Intersection Point along segment:
            % DS    = (XRAY-X2_1) / (X2_2-X2_1)
            % SVARI = X3_1 + DS*(X3_2-X3_1)
            SVARI = X3_1 + (XRAY-X2_1) * (X3_2-X3_1) / (X2_2-X2_1);
            % LOW and HIGH media type, using the segment definition:
            ICRSI(LOW_IND) = GAM(LOW_IND)*ISSEG(LOW_IND) + GAM(HIGH_IND)*ISSEG(HIGH_IND);
            ICRSI(HIGH_IND)= GAM(LOW_IND)*ISSEG(HIGH_IND)+ GAM(HIGH_IND)*ISSEG(LOW_IND);
            ICRSI(HIGH_IND+1)=IBM_UNDEFINED;
            SCRSI = ISEG;
            % Insertion sort:
            [ierr2]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI);
            continue
         end
         disp(['Error GET_IS_SOLID_PT NVEC(JAXIS): Missed segment=' num2str(ISEG)])
      end
    end 
         
end

% Do we have any intersections?
if ( IBM_N_CRS == 0 ); IS_SOLID =false; return; end

[ierr]=COLLAPSE_CROSSINGS(BODINT_PLANE2,X1AXIS,X2AXIS,X3AXIS,XRAY,X1PLN,2);

[IS_GASPHASE]=GET_IS_GASPHASE(SCEN);

IS_SOLID = ~IS_GASPHASE;

return