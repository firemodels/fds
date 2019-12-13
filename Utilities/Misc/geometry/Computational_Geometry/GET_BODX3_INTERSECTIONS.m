function [ierr]=GET_BODX3_INTERSECTIONS(X2AXIS,X3AXIS,X2LO,X2HI)

global GEOMEPS NOD1 NOD2 IAXIS JAXIS KAXIS BODINT_PLANE
global X2FACE DX2FACE
global X1NOC X2NOC X3NOC

ierr=1;

for ISEG=1:BODINT_PLANE.NSEGS

   if (BODINT_PLANE.X3ALIGNED(ISEG)); continue; end % This segment is not aligned with x3.

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
   STANI(IAXIS:JAXIS) = [ (X2_2-X2_1), (X3_2-X3_1) ]*SLEN^(-1.);

   MINX = min(X2_1,X2_2);
   MAXX = max(X2_1,X2_2);
   if (X2NOC==0)
      % Optimized for UG:
      JSTR = max(X2LO, ceil((  MINX-GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO))+X2LO);
      JEND = min(X2HI,floor((  MAXX+GEOMEPS-X2FACE(X2LO))/DX2FACE(X2LO))+X2LO);
   else
      if ((MINX-GEOMEPS-X2FACE(X2LO)) < 0.)
         JSTR=X2LO;
      elseif((MINX-GEOMEPS-X2FACE(X2HI)) >= 0.)
         JSTR=X2HI+1;
      else
         for JJ=X2LO:X2HI
            if ((MINX-GEOMEPS-X2FACE(JJ)) >= 0. && (MINX-GEOMEPS-X2FACE(JJ+1)) < 0. )
               JSTR = JJ+1;
               break
            end
         end
      end
      if ((MAXX+GEOMEPS-X2FACE(X2LO)) < 0.)
         JEND=X2LO-1;
      elseif((MAXX+GEOMEPS-X2FACE(X2HI)) >= 0.)
         JEND=X2HI;
      else
         for JJ=X2LO:X2HI
            if ((MAXX+GEOMEPS-X2FACE(JJ)) >= 0. && (MAXX+GEOMEPS-X2FACE(JJ+1)) < 0. )
               JEND = JJ;
               break
            end
         end
      end
   end

   for JJ=JSTR:JEND

      % S coordinate along segment:
      DX2_1 = X2_2 - X2FACE(JJ);
      DX2_2 = X2FACE(JJ) - X2_1;
      XI1   = DX2_1 / (X2_2-X2_1);
      XI2   = DX2_2 / (X2_2-X2_1);
      DV(IAXIS:JAXIS) = [ DX2_2, (XI1-1.)*X3_1+XI2*X3_2 ];
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

end

ierr=0;

return