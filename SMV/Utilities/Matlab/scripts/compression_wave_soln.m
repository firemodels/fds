% McDermott
% 4-8-10
% compression_wave_soln.m

function [rho] = compression_wave_soln(rho0,x,y,a,c,t);

b = sqrt(-1+a^2);
d = sqrt(-1+c^2);

x0 = 2*atan( b/a*tan( atan( (1+a*tan(x/2))/b ) - b*t/2 ) - 1/a );
y0 = 2*atan( d/c*tan( atan( (1+c*tan(y/2))/d ) - d*t/2 ) - 1/c );

Ix0 = log( -a^2 - cos( 2*atan( (1+a*tan(x0/2))/b ) ) + b*sin( 2*atan( (1+a*tan(x0/2))/b ) ) );
Iy0 = log( -c^2 - cos( 2*atan( (1+c*tan(y0/2))/d ) ) + d*sin( 2*atan( (1+c*tan(y0/2))/d ) ) );

Ix = log( -a^2 - cos(b*t + 2*atan( (1+a*tan(x0/2))/b ) ) + b*sin(b*t + 2*atan( (1+a*tan(x0/2))/b ) ) );
Iy = log( -c^2 - cos(d*t + 2*atan( (1+c*tan(y0/2))/d ) ) + d*sin(d*t + 2*atan( (1+c*tan(y0/2))/d ) ) );

q0 = log(rho0);
q = q0 + Ix-Ix0 + Iy-Iy0;

rho = exp(q);