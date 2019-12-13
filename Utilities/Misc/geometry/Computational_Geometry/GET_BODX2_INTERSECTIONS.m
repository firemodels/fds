function [ierr]=GET_BODX2_INTERSECTIONS(X2AXIS,X3AXIS,X3RAY)

global GEOMEPS NOD1 NOD2 IAXIS JAXIS KAXIS BODINT_PLANE
global IBM_N_CRS IBM_SEG_CRS IBM_SVAR_CRS IBM_SEG_TAN


ierr=1;

if ( BODINT_PLANE.NSEGS == 0); return; end

for ICRS=1:IBM_N_CRS

   ISEG = IBM_SEG_CRS(ICRS);

   if (ISEG < 0); continue; end % it is a single point element.

   SEG(NOD1:NOD2)    = BODINT_PLANE.SEGS(NOD1:NOD2,ISEG);
   XYZ1(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD1));
   XYZ2(IAXIS:KAXIS) = BODINT_PLANE.XYZ(IAXIS:KAXIS,SEG(NOD2));

   % x2_x3 of segment point 1:
   X2_1 = XYZ1(X2AXIS); X3_1 = XYZ1(X3AXIS);
   % x2_x3 of segment point 2:
   X2_2 = XYZ2(X2AXIS); X3_2 = XYZ2(X3AXIS);

   % Segment length:
   SLEN = sqrt( (X2_2-X2_1)^2. + (X3_2-X3_1)^2. );

   % Unit vector along segment:
   STANI(IAXIS:JAXIS) = IBM_SEG_TAN(IAXIS:JAXIS,ICRS);

   % S coordinate along segment:
   DV(IAXIS:JAXIS)= [(IBM_SVAR_CRS(ICRS)-X2_1), (X3RAY-X3_1)];
   SBOD = DV(IAXIS)*STANI(IAXIS)+DV(JAXIS)*STANI(JAXIS);

   % If crossing is at one end, cycle:
   if ( (abs(SBOD) < GEOMEPS) || (abs(SBOD-SLEN) < GEOMEPS) ); continue; end

   % Add crossing to BODINT_PLANE, insertion sort:
   NBCROSS  = BODINT_PLANE.NBCROSS(ISEG) + 1;
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