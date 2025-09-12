function u_new = project_cells(u, u1, dt, FA_cell, FB_cell, gradH, stage)
% FA_cell     : Nc x 2  (advective + body) cell vector
% divTau_cell : Nc x 2  (∇·τ) at cells
% rho         : Nc x 1
% FB_cell     : Nc x 2  (-p ∇(1/ρ)) at cells
% gradH       : Nc x 2  ∇H at cells
% Output:
%   u_new     : Nc x 2
if(stage==1)
   u_new = u - dt*( FA_cell + FB_cell + gradH );
else
   u_new = 0.5*(u+u1) - dt/2*( FA_cell + FB_cell + gradH );
end
