function T = slab(L,x,k,rho,cp,h,t,T0,Tinf,N,options)
% SLAB    Exact solution of the temperature inside a slab
%   with initial temperature T0 and both surfaces 
%   subject to heat transfer from Tinf temperature from t = 0.
%
% T = slab(L,x,k,rho,cp,h,t,T0,Tinf,N)
%
%   L       Slab half-thickness (one side) (m)
%   x       distance from the slab center (m)
%   k       conductivity    (W/Km)
%   rho     density         (kg/m3)
%   cp      heat capasity   (J/kgK)
%   h       heat transfer coefficient (W/m2K)
%   t       time vector             (s)
%   T0      Initial temperature     (C)
%   Tinf    Gas temperature         (C)
%   N       Number of terms in the series
%   options (optional) options for fzero

if (nargin<11)
    options = optimset('fzero');
end
Bi = h*L/k;
alpha = k/(rho*cp);
% roots
for n = 1:N
    la(n) = fzero(@(x) cot(x*L)-x*L/Bi,[n-1+0.00001 n-0.00001]*pi/L,options);
end
%In Drysdale, it was
%for n = 1:N
%    la(n) = fzero(@(x) cot(x)-x*L/Bi,[n-1+0.001 n-0.001]*pi,options);
%end
% test roots
%max(cot(la*L)-la*L/Bi)
%
laL = la*L;
for n = 1:N
    laN = laL(n);
    Sla(n) = sin(laN)/(laN + sin(laN)*cos(laN)) ;
end
%
theta0 = T0-Tinf;
R(length(t),length(x)) = 0;
for i = 1:length(t)
    if (t(i)==0)
        R(1,:) = ones([1 length(x)]);
    else
    for n = 1:N
        S = Sla(n) * exp(-la(n)^2*alpha*t(i));
        for j = 1:length(x)
            R(i,j) = R(i,j) + 2 * S * cos(la(n)*x(j));
        end
    end
    end
end
theta = R*theta0;
T = theta + Tinf;
 