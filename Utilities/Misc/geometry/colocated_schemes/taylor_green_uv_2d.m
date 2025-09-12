function UV = taylor_green_uv_2d(x, y, t, U0, kx, ky, nu)

global ux0 uy0

% x,y can be vectors or arrays (same size). Returns [u,v] stacked along 2nd dim.
k2 = kx*kx + ky*ky;
decay = exp(-nu*k2*t);
u = -U0 * cos(kx*(x-ux0*t)) .* sin(ky*(y-uy0*t)) * decay + ux0;
v =  U0 * sin(kx*(x-ux0*t)) .* cos(ky*(y-uy0*t)) * decay + uy0;
UV = [u(:), v(:)];
end