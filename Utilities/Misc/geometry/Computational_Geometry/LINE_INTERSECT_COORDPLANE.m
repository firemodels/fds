function [XYZ_INT,INTFLG]=LINE_INTERSECT_COORDPLANE(X1AXIS,X1PLN,PLNORMAL,LNC)

global IAXIS JAXIS KAXIS NOD1 NOD2 GEOMEPS

% Initialize:
INTFLG = false;
XYZ_INT(IAXIS:KAXIS) = 0.;

% Preliminary calculations:
DVEC(IAXIS:KAXIS) = LNC(IAXIS:KAXIS,NOD2) - LNC(IAXIS:KAXIS,NOD1);
NMDV  = sqrt( DVEC(IAXIS)^2. + DVEC(JAXIS)^2. + DVEC(KAXIS)^2. );
DIRV  = DVEC(IAXIS:KAXIS) * NMDV^(-1.);
DENOM = DIRV(IAXIS)*PLNORMAL(IAXIS) +DIRV(JAXIS)*PLNORMAL(JAXIS) +DIRV(KAXIS)*PLNORMAL(KAXIS);
PLNEQ = LNC(IAXIS,NOD1)*PLNORMAL(IAXIS) + ...
        LNC(JAXIS,NOD1)*PLNORMAL(JAXIS) + ...
        LNC(KAXIS,NOD1)*PLNORMAL(KAXIS) - X1PLN;

% Line parallel to plane:
if ( abs(DENOM) < GEOMEPS )
   % Check if seg lies on plane or not.
   % Do this by checking if node one of segment is on plane.
   if ( abs(PLNEQ) < GEOMEPS )
      XYZ_INT(IAXIS:KAXIS) = LNC(IAXIS:KAXIS,NOD1); XYZ_INT(X1AXIS) = X1PLN;
      INTFLG = true;
   end
   return
end

% Non parallel case:
TLINE = -PLNEQ/DENOM;  % Coordinate along the line LNC.
XYZ_INT(IAXIS:KAXIS,1) = LNC(IAXIS:KAXIS,NOD1)' + TLINE*DIRV(IAXIS:KAXIS); % Intersection point.
XYZ_INT(X1AXIS,1) = X1PLN; % Force X1AXIS coordinate to be the planes value.
INTFLG = true;

return
