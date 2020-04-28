function [IS_GASPHASE]=GET_IS_GASPHASE(SCEN)

global IBM_N_CRS IBM_SVAR_CRS IBM_IS_CRS GEOMEPS IBM_GS IBM_SG

% Count GS,SG intersections from both sides:
IS_GASPHASE_LEFT = true;
for ICRS=1:IBM_N_CRS
   if (SCEN < IBM_SVAR_CRS(ICRS)-GEOMEPS/2.); continue; end

   % If solid change state:
   if ( (IBM_IS_CRS(ICRS) == IBM_GS) || (IBM_IS_CRS(ICRS) == IBM_SG) )
        IS_GASPHASE_LEFT = ~IS_GASPHASE_LEFT;
   end
end

IS_GASPHASE_RIGHT = true;
for ICRS=IBM_N_CRS:-1:1
    if (SCEN > IBM_SVAR_CRS(ICRS)+GEOMEPS/2.); continue; end

    % If solid change state:
    if ( (IBM_IS_CRS(ICRS) == IBM_GS) || (IBM_IS_CRS(ICRS) == IBM_SG) )
        IS_GASPHASE_RIGHT = ~IS_GASPHASE_RIGHT;
    end
end

% If at least one of left and right are true -> add
% IBM_GASPHASE cut-edge:
IS_GASPHASE = IS_GASPHASE_LEFT || IS_GASPHASE_RIGHT;

return