function [IND_PI,BODINT_PLANE]=GET_BODINT_NODE_INDEX(BODINT_PLANE,X2AXIS,X3AXIS,XYZ)

global IAXIS KAXIS GEOMEPS

% Test if XYZ is already on BODINT_PLANE.XYZ:
IND_PI = -1; % Initialize to negative index.
% Linear Search:
inlist=false;
for INOD=1:BODINT_PLANE.NNODS
  DIFFX2 = BODINT_PLANE.XYZ(X2AXIS,BODINT_PLANE.NOD_PERM(INOD))-XYZ(X2AXIS);
  if ( DIFFX2 > GEOMEPS )
     break
  elseif ( abs(DIFFX2) <= GEOMEPS)
     DIFFX3 = BODINT_PLANE.XYZ(X3AXIS,BODINT_PLANE.NOD_PERM(INOD))-XYZ(X3AXIS);
     if ( DIFFX3 > GEOMEPS )
        break
     elseif ( abs(DIFFX3) <= GEOMEPS )
        IND_PI = BODINT_PLANE.NOD_PERM(INOD);
        inlist = true;
        return
     end
  end
end

   
if(~inlist)
    
    if (BODINT_PLANE.NNODS == 0); INOD=1; end
    
    % Insert add NOD_PERM permutation array, O(NP) operation:
    for INOD2=BODINT_PLANE.NNODS+1:-1:INOD+1
        BODINT_PLANE.NOD_PERM(INOD2) = BODINT_PLANE.NOD_PERM(INOD2-1);
    end
    
    IND_PI = BODINT_PLANE.NNODS + 1;
    BODINT_PLANE.NNODS = IND_PI;
    BODINT_PLANE.NOD_PERM(INOD) = IND_PI;
    BODINT_PLANE.XYZ(IAXIS:KAXIS,IND_PI) = XYZ(IAXIS:KAXIS);
    
end

return