% Antonellis
% 06-23-2010
% section2_soln.m

function [rho] = section2_soln(rho0,x,y,B,w,t);
   
   x0 = 2*atan(tan(x/2)*exp(-B/w*sin(w*t)));
   y0 = 2*atan(tan(y/2)*exp(-B/w*sin(w*t)));

   q0 = log(rho0);
   q = q0 + log((1 + (tan(x0/2)).^2.*exp(2*B/w*sin(w*t)))./(1+(tan(x0/2)).^2)) ...
          + log((1 + (tan(y0/2)).^2.*exp(2*B/w*sin(w*t)))./(1+(tan(y0/2)).^2)) ...
          - 2*B/w*sin(w*t); % q(x,y,t) in Verification Guide
                                     
   rho = exp(q);