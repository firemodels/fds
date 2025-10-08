function dr = periodic_delta_2d(dx_raw, Lx, Ly)
% Wrap dx_raw = [dx, dy] into (-L/2, L/2]
dr = dx_raw;
% x
if     dr(1) >  0.5*Lx, dr(1) = dr(1) - Lx;
elseif dr(1) <= -0.5*Lx, dr(1) = dr(1) + Lx; end
% y
if     dr(2) >  0.5*Ly, dr(2) = dr(2) - Ly;
elseif dr(2) <= -0.5*Ly, dr(2) = dr(2) + Ly; end
end