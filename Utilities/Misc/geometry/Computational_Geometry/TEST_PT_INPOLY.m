function [PTSFLAG]=TEST_PT_INPOLY(NP,XY,XY1)

global IAXIS JAXIS GEOMEPS

PTSFLAG = false;
RCROSS  = 0;
LCROSS  = 0;

% ADD first point location at the end of XY (assumes IBM_MAXVERTS_FACE > NP):
XY(IAXIS:JAXIS,NP+1) = XY(IAXIS:JAXIS,1);

% Shift origin to XY1:
for IP=1:NP+1
    XY(IAXIS:JAXIS,IP) = XY(IAXIS:JAXIS,IP) - XY1(IAXIS:JAXIS)';
end

% For each edge test against rays x=0, y=0:
for IP=1:NP
   % Check if edges first point is vertex:
   if ( (abs(XY(IAXIS,IP)) < GEOMEPS) && ...
        (abs(XY(JAXIS,IP)) < GEOMEPS) )
      PTSFLAG = true;
      return
   end
   % Check if edge crosses x axis:
   RS = (XY(JAXIS,IP) > 0.) ~= (XY(JAXIS,IP+1) > 0.);
   LS = (XY(JAXIS,IP) < 0.) ~= (XY(JAXIS,IP+1) < 0.);

   if ( RS || LS )
      % Intersection:
      XPT =(XY(IAXIS,IP  )*XY(JAXIS,IP+1) - ...
            XY(JAXIS,IP  )*XY(IAXIS,IP+1)) / (XY(JAXIS,IP+1)-XY(JAXIS,IP));

      if (RS && (XPT > 0.)); RCROSS = RCROSS + 1; end
      if (LS && (XPT < 0.)); LCROSS = LCROSS + 1; end
   end
end

if ( mod(RCROSS,2) ~= mod(LCROSS,2) ) % Point on edge
   PTSFLAG = true;
   return
end

if ( mod(RCROSS,2) == 1) % Point inside
   PTSFLAG = true;
   return
end

return
