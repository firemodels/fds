% PARK HYUN WOOK, Yonsei University
% 8-15-2012
% blasius_analytic.m

function [eta,fp] = blasius_analytic(u0, zmax)

% etamax = maximum eta for calculation
% steps = number of steps between 0 and etamax
% fppwall = initial (wall) value of 2nd derivative of Blasius function
% outputs:
% eta - the similarity coordinate normal to the wall
% f, fp, fpp, fppp - the Blasius function and it first 3 derivatives

etamax=zmax/sqrt(1E-3/1.199*0.05/u0);
steps=257;
deta = etamax/(steps-1);
eta = zeros(steps,1);
f = zeros(steps,1);
fp = zeros(steps,1);
fpp = zeros(steps,1);
fppp = zeros(steps,1);
% initial guess for fpp
fpp(1) = 0.3318;
for i=1:steps-1
    eta(i+1) = eta(i) + deta;
    % predictor
    %1st
    k1(1)=fp(i);
    k1(2)=fpp(i);
    k1(3)=-f(i)*fpp(i)/2;
    
    fbar = f(i) + 0.5*deta * k1(1);
    fpbar = fp(i) + 0.5*deta * k1(2);
    fppbar = fpp(i) + 0.5*deta * k1(3);
    
    %2nd
    k2(1)=fpbar;
    k2(2)=fppbar;
    k2(3)=-fbar*fppbar/2;
    
    fbar = f(i) + 0.5*deta * k2(1);
    fpbar = fp(i) + 0.5*deta * k2(2);
    fppbar = fpp(i) + 0.5*deta * k2(3);
    
    %3rd
    k3(1)=fpbar;
    k3(2)=fppbar;
    k3(3)=-fbar*fppbar/2;
    
    fbar = f(i) + deta * k3(1);
    fpbar = fp(i) + deta * k3(2);
    fppbar = fpp(i) + deta * k3(3);
    
    %4th
    k4(1)=fpbar;
    k4(2)=fppbar;
    k4(3)=-fbar*fppbar/2;
    
    % corrector
    f(i+1) = f(i) + deta * (k1(1)+2*k2(1)+2*k3(1)+k4(1))/6;
    fp(i+1) = fp(i) + deta * (k1(2)+2*k2(2)+2*k3(2)+k4(2))/6;
    fpp(i+1) = fpp(i) + deta * (k1(3)+2*k2(3)+2*k3(3)+k4(3))/6;
    fppp(i+1) = -f(i+1)*fpp(i+1)/2;
end
