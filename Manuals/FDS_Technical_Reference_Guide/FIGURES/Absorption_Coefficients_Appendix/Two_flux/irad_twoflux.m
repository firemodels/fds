function [q,X] = irad_twoflux(L,qin,kappa,T,N)
% [q,X] = irad_twoflux(L,qin,kappa,T[,N]);
% 
% Function irad_twoflux returns the 
% transmitted radiative heat flux q at depths X, using two-flux method.
% Inputs:   L      layer depth (m)
%           qin    incident heat flux (kW/m2)
%           kappa  appropriate mean absorption coefficient (1/m)
%           T      material temperature (C)
%           N      number of grid cells (optional)
if (nargin<5), N = 20; end
sigma = 5.67E-8;
T = T+273.15;
Ib = sigma*T^4;
%
% linear grid
dx = L/N;
X = [1:N]*dx;
%
dx = X(1);
q(1) = (qin + 2*dx*kappa*Ib)/(1+2*dx*kappa);
for i = 2:N
    dx = X(i)-X(i-1);
    q(i) = (q(i-1) + 2*dx*kappa*Ib)/(1+2*dx*kappa); 
end
