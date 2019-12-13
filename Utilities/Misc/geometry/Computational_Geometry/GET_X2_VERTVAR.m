function [ierr]=GET_X2_VERTVAR(X1AXIS,X2LO,X2HI,NM,I,KK)

global IBM_GG IBM_SS IBM_GS IBM_SG IAXIS GEOMEPS IBM_SOLID IBM_IS_CRS IBM_SVAR_CRS MESHES
global X2FACE DX2FACE
global IBM_VGSC
global IBM_N_CRS

ierr=1;

X2NOC=0;

% Work By Edge, Only one X1AXIS=IAXIS needs to be used:
switch(X1AXIS)
case(IAXIS)
   X2LO_LOC = X2LO;
   X2HI_LOC = X2HI;
   % Case of GG, SS points:
   for ICRS=1:IBM_N_CRS
      % If is_crs(icrs) == GG, SS, SGG see if crossing is
      % exactly on a Cartesian cell vertex:
      switch(IBM_IS_CRS(ICRS))
      case{IBM_GG, IBM_SS}
         JSTR = X2LO_LOC; JEND = X2HI_LOC;
         if(X2NOC==0) 
            % Optimized and will ONLY work for Uniform Grids:
            JSTR = max(X2LO_LOC, floor((IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(X2LO_LOC))/DX2FACE(X2LO_LOC)) + X2LO_LOC);
            JEND = min(X2HI_LOC, ceil((IBM_SVAR_CRS(ICRS)+GEOMEPS-X2FACE(X2LO_LOC))/DX2FACE(X2LO_LOC)) + X2LO_LOC);
         end

         for JJ=JSTR:JEND
            % Crossing on Vertex?
            if ( abs(X2FACE(JJ)-IBM_SVAR_CRS(ICRS)) < GEOMEPS ) 
               MESHES(NM).VERTVAR(I,JJ,KK,IBM_VGSC) = IBM_SOLID;
               break
            end
         end

      end
   end

   % Other cases:
   for ICRS=1:IBM_N_CRS-1
      % Case GS-SG: All Cartesian vertices are set to IBM_SOLID.
      if (IBM_IS_CRS(ICRS) == IBM_GS) 
         % Find corresponding SG intersection:
         for ICRS1=ICRS+1:IBM_N_CRS
            if (IBM_IS_CRS(ICRS1) == IBM_SG); break; end
         end
         JSTR = X2LO_LOC; JEND = X2HI_LOC;
         if(X2NOC==0) 
            % Optimized for UG:
            JSTR = max(X2LO_LOC, ceil(( IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(X2LO_LOC))/DX2FACE(X2LO_LOC)) + X2LO_LOC);
            JEND = min(X2HI_LOC, floor((IBM_SVAR_CRS(ICRS1)+GEOMEPS-X2FACE(X2LO_LOC))/DX2FACE(X2LO_LOC)) + X2LO_LOC);
         else
            if ((IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(X2LO_LOC))    <  0.) 
               JSTR=X2LO_LOC;
            elseif((IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(X2HI_LOC)) >= 0.) 
               JSTR=X2HI_LOC+1;
            else
               for JJ=X2LO_LOC:X2HI_LOC
                  if((IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(JJ))   >= 0. && ...
                     (IBM_SVAR_CRS(ICRS)-GEOMEPS-X2FACE(JJ+1)) <  0. ) 
                     JSTR = JJ+1;
                     break
                  end
               end
            end
            if ((IBM_SVAR_CRS(ICRS1)+GEOMEPS-X2FACE(X2LO_LOC)) < 0.) 
               JEND=X2LO_LOC-1;
            elseif((IBM_SVAR_CRS(ICRS1)+GEOMEPS-X2FACE(X2HI))  >= 0.) 
               JEND=X2HI_LOC;
            else
               for JJ=X2LO_LOC:X2HI_LOC
                  if((IBM_SVAR_CRS(ICRS1)+GEOMEPS-X2FACE(JJ))   >= 0. && ...
                     (IBM_SVAR_CRS(ICRS1)+GEOMEPS-X2FACE(JJ+1)) <  0. ) 
                     JEND = JJ;
                     break
                  end
               end
            end
         end

         for JJ=JSTR:JEND
            MESHES(NM).VERTVAR(I,JJ,KK,IBM_VGSC) = IBM_SOLID;
         end
      end
   end
end

ierr = 0;

return