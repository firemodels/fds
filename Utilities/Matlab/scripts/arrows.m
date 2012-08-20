%-------------------------------------------------------------------------%
% FILE: arrows.m
% AUTHOR: Max Gould
% LAST DATE EDITED: 27 June 2012
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Generate filled arrow heads for quiver plot.
% Should specify ('ShowArrowHead','off') in quiver argument.
%
% [] = arrows(x,y,u,v,[alpha],[beta])
%
% Arguments:
%     x - array of x-coordinates for every point in the vector field
%
%     y - array of y-coordinates for every point in the vector field
%
%     u - array of x-component vector magnitudes
%         for every point in the vector field
%
%     v - array of y-component vector magnitudes
%         for every point in the vector field
%
%    [optional] alpha - specify primary scaling factor,
%                       determines relative size of arrow heads
%
%    [optional] beta - specify secondary scaling factor,
%                      determines width to length ratio of arrow heads
%-------------------------------------------------------------------------%

function [] = arrows(x,y,u,v,varargin)

%scales arrow-head placement with vector lines
scale = 1.6*abs(x(1,1)-x(1,2));

% Convert arrays to vectors and scale
x=x(:).';
y=y(:).';
u=scale*u(:).';
v=scale*v(:).';

% Check array dimensions
if length(x) ~= length(y)
    error('myApp:argChk',['Plot range dimensions do not match. \n',...
          'Try using [x,y]=meshgrid'])
end
if length(u) ~= length(v)
    error('myApp:argChk','Vector array dimensions do not match')
end

% Define arrow-head dimensions
    % alpha = primary scaling factor (arrow size)
    % beta  = secondary scaling factor (width to length ratio)
if nargin == 4
    % Define default values
    alpha = 0.32;
    beta  = 0.25;
elseif nargin == 5
    alpha = varargin{1};
    beta  = 0.25;
elseif nargin ==6
    alpha = varargin{1};
    beta = varargin{2};
elseif nargin > 6
    error('myApp:argChk','Unexpected number of arguments')
end

% Allocate space for arrays
hu = zeros(length(u),4);
hv = zeros(length(u),4);

% Define arrow geometry and plot arrows
for i = 1:length(u)
    hu(i,1) = x(i)+u(i)-alpha*(u(i)+beta*(v(i)+eps));
    hu(i,2) = x(i)+u(i);
    hu(i,3) = x(i)+u(i)-alpha*(u(i)-beta*(v(i)+eps));
    hu(i,4) = x(i)+u(i)-alpha*(u(i)+beta*(v(i)+eps));
    hv(i,1) = y(i)+v(i)-alpha*(v(i)-beta*(u(i)+eps));
    hv(i,2) = y(i)+v(i);
    hv(i,3) = y(i)+v(i)-alpha*(v(i)+beta*(u(i)+eps));
    hv(i,4) = y(i)+v(i)-alpha*(v(i)-beta*(u(i)+eps));
    fill(hu(i,:),hv(i,:),[0.3 0.7 1]);
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%