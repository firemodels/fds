function [up, vp] = taylor_green_uv_2d_point(xpi, ypi, t, U0, kx, ky, nu, ux0, uy0)
% TAYLOR_GREEN_UV_2D_POINT  TG velocities at a single (xp,yp).
if nargin < 9, uy0 = 0; end
if nargin < 8, ux0 = 0; end
xp = xpi-ux0*t;
yp = ypi-uy0*t;
k2    = kx*kx + ky*ky;
decay = exp(-nu * k2 * t);
up = -U0 * cos(kx*xp) * sin(ky*yp) * decay + ux0;
vp =  U0 * sin(kx*xp) * cos(ky*yp) * decay + uy0;
end
