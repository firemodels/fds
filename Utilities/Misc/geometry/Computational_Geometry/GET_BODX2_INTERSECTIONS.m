function [ierr]=GET_BODX2_INTERSECTIONS(X2AXIS,X3AXIS,X3RAY)

global GEOMEPS NOD1 NOD2 IAXIS JAXIS KAXIS BODINT_PLANE
global IBM_N_CRS IBM_SEG_CRS IBM_SVAR_CRS IBM_SEG_TAN


ierr=1;

if ( BODINT_PLANE.NSEGS == 0); return; end


for ISEG=1:BODINT_PLANE.NSEGS

   if (BODINT_PLANE.X2ALIGNED(ISEG)); continue; end % This segment is aligned with x2.

   SEG(NOD1:NOD2)    = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD1));
   XYZ2(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD2));

   % x2_x3 of segment point 1:
   X2_1 = XYZ1(X2AXIS); X3_1 = XYZ1(X3AXIS);
   % x2_x3 of segment point 2:
   X2_2 = XYZ2(X2AXIS); X3_2 = XYZ2(X3AXIS);

   if( (X3RAY-max(X3_1,X3_2)) > GEOMEPS); continue; end
   if( (min(X3_1,X3_2)-X3RAY) > GEOMEPS); continue; end
   
   % Segment length:
   SLEN = sqrt( (X2_2-X2_1)^2. + (X3_2-X3_1)^2. );

   % Unit vector along segment:
   STANI(IAXIS:JAXIS) = [ (X2_2-X2_1), (X3_2-X3_1) ]*SLEN^(-1.);

   % S coordinate along segment:
   DX3_1 = X3_2 - X3RAY;
   DX3_2 = X3RAY- X3_1;
   XI1   = DX3_1 / (X3_2-X3_1);
   XI2   = DX3_2 / (X3_2-X3_1);
   DV(IAXIS:JAXIS) = [(XI1-1.)*X2_1+XI2*X2_2, DX3_2];
   SBOD = DV(IAXIS)*STANI(IAXIS)+DV(JAXIS)*STANI(JAXIS);

   % If crossing is already defined, cycle:
   NBCROSS = BODINT_PLANE.NBCROSS(ISEG);
   ISCONT = false;
   for IBCR=1:NBCROSS
       if ( abs(SBOD-BODINT_PLANE.SVAR(IBCR,ISEG)) < GEOMEPS )
           ISCONT = true;
           break
       end
   end
   if (ISCONT); continue; end
   
   % Add crossing to BODINT_PLANE, insertion sort:
   NBCROSS = BODINT_PLANE.NBCROSS(ISEG) + 1;
   
   % Test-reallocate BODINT_PLANE.SVAR
   BODINT_PLANE.SVAR(NBCROSS,ISEG) = 1./GEOMEPS;
   for IBCR=1:NBCROSS
       if ( SBOD < BODINT_PLANE.SVAR(IBCR,ISEG) ); break; end
   end
   
   % Here copy from the back (updated nbcross) to the ibcr location:
   for IDUM = NBCROSS:-1:IBCR+1
       BODINT_PLANE.SVAR(IDUM,ISEG) = BODINT_PLANE.SVAR(IDUM-1,ISEG);
   end
   BODINT_PLANE.SVAR(IBCR,ISEG) = SBOD;
   BODINT_PLANE.NBCROSS(ISEG)   = NBCROSS;

end

ierr=0;

return