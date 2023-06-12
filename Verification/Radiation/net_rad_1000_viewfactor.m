function phi = viewfactor(x,b)
%a = sqrt(2)*b;
alfa = pi/4;
phi=0.5*(1+(b*cos(alfa)-x)./(x.^2-2*b*x*cos(alfa)+b^2).^(1/2));
end