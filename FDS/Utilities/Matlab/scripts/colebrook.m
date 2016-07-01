% McDermott
% 2-23-09
% colebrook.m
%
% input:
% Re = Reynolds number
% RR = relative roughness
% ff = a guess for the value of f
% tol = convergence tolerance (relative error between f and ff)
%
% output:
% f = Moody friction factor based on Colebrook formula
%
% ref:
% Munson, Young and Okiishi, Fundamentals of Fluid Mechanics, Wiley, 1990.

function [f,error,iter] = colebrook(Re,RR,ff,tol)
iter=0;
error=tol+1;
while error > tol
    iter=iter+1;
    f = (-2*log10(RR/3.7 + 2.51/(Re*sqrt(ff))))^(-2);
    error = abs((f-ff)/ff);
    ff = f; % Picard iteration
end
