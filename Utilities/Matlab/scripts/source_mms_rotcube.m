% Immersed rotated square object manufactured solution:
% 
% Velocity fields and pressure in local axes are based on Weinan and Liu 
% (CMS 2003) test case:
% Marcos Vanella - 2019.
%% ------------------------------------------------------------------------
close all
clear all
clc

syms x y k nu t real

%% Velocities and pressure manufactured solution forcing term:
% Cartesian Case:
u  = -sin(t)*sin(k*x)^2*sin(2*k*y);
v  =  sin(t)*sin(2*k*x)*sin(k*y)^2;
p  =  sin(t)/4*(2+cos(2*k*x))*(2+cos(2*k*y))-sin(t);


% Diffusive and advective terms Fv = - nu Lap(u) + u.Grad(u):
Fvx = simplify(-nu*(diff(u,x,2)+diff(u,y,2)) + diff(u*u,x) + diff(u*v,y));
Fvy = simplify(-nu*(diff(v,x,2)+diff(v,y,2)) + diff(v*u,x) + diff(v*v,y));

dpdx = simplify(diff(p,x));
dpdy = simplify(diff(p,y));

dudt   = simplify(diff(u,t));
dvdt   = simplify(diff(v,t));


% Forcing term: Udot = -(Fv+Grad(p)+fe_v) -> fe_v = -(Fv+Grad(p))-Udot
fex_v = simplify(-(Fvx+dpdx)-dudt);
fey_v = simplify(-(Fvy+dpdy)-dvdt);

% Divergence of U:
divU = simplify(diff(u,x) + diff(v,y));


%% Scalar function manufactured solution:
syms rho D_Z A_Z Z_mean gam real

% Manufactured solution for scalar Z(x,y,t):
Z=A_Z/3*sin(t)*((1-cos(2*k*(x-gam)))*(1-cos(2*k*(y-gam)))-1)+Z_mean;

% d(rho*Z)/dt:
drhoZdt = simplify(diff(rho*Z,t));

% Diffusive fluxes for scalar:
diffZx = simplify(-rho*D_Z*diff(Z,x));
diffZy = simplify(-rho*D_Z*diff(Z,y));

% Advective and diffusive terms for scalar:
Fz = simplify(diff(rho*Z*u,x)+diff(rho*Z*v,y)+ ...
              diff(diffZx,x)+diff(diffZy,y));
           
% Forcing for Z eqn.: d(rho*Z)/dt=-(Fz+fe_Z) -> fe_Z = -F_z-d(rho*Z)/dt; 
fe_Z=-simplify(drhoZdt + Fz);

return 

% Plot functions in Eulerian reference frame:
IAXIS=1; JAXIS=2;
% Scalar Function Parameters:
k        =   1;       % Wave number.
A_Z      = 0.1;       % Scalar function amplitude.
Z_mean   = 0.15;      % scalar function mean.
rotangle = atan(1/2); % Rotation angle can be atan(1/r) where r=1,2,3,4,...
gam      = pi/2;      % Shift for scalar function.

% Domain size:
DELY=pi;
if(abs(rotangle) < 1.e-12)
    DELX=0;
else
    DELX=DELY/tan(rotangle);
end
Len=sqrt(DELX^2+DELY^2);

xmin=pi-Len;
xmax=pi+Len;
ymin=pi-Len;
ymax=pi+Len;

Lx=xmax-xmin;
Ly=ymax-ymin;

% Define arrays:
ncell=192; af=2;
nface=ncell+1;
xv=linspace(xmin,xmax,nface);
xc=0.5*(xv(1:ncell)+xv(2:nface));
yv=linspace(ymin,ymax,nface);
yc=0.5*(yv(1:ncell)+yv(2:nface));

u  =zeros(nface,ncell);
v  =zeros(ncell,nface);
P  =zeros(ncell,ncell); Z=P;
div=zeros(ncell,ncell);
uc =zeros(ncell,ncell);
vc =zeros(ncell,ncell);

displ=pi*[1. 1.]';
dispb=pi*[1. 1.]';
dispx=pi/2; dispy=dispx;

RotMatrix=[cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];

nt=25;
tv=linspace(0,2*pi,nt);

% Body for plotting:
xbdl=[0 pi pi  0 0];
ybdl=[0  0 pi pi 0];

for ip=1:length(xbdl)    
   xbdg   = RotMatrix*([xbdl(ip) ybdl(ip)]'-[dispx dispy]') + dispb;
   xbd(ip)= xbdg(IAXIS);
   ybd(ip)= xbdg(JAXIS);
end

% figure:
figure
hold on
it=0;
for t=tv 
    it=it+1;
    
    % U velocity:
    for j=1:ncell
        for i=1:nface      
            xglob=[xv(i) yc(j)]';
            xloc=RotMatrix'*(xglob-displ)+[dispx dispy]';
            if ( (xloc(IAXIS) > 0) && (xloc(IAXIS) < pi) && ...
                 (xloc(JAXIS) > 0) && (xloc(JAXIS) < pi) )
             u(i,j)= 0;
            else            
             uloc(IAXIS)=-sin(t)*sin(k*xloc(IAXIS))^2*sin(2*k*xloc(JAXIS));
             uloc(JAXIS)= sin(t)*sin(2*k*xloc(IAXIS))*sin(k*xloc(JAXIS))^2;
             uglob= RotMatrix * uloc';
             u(i,j)= uglob(IAXIS);
            end
        end
    end    
   
    % V velocity:
    for j=1:nface
        for i=1:ncell           
            xglob=[xc(i) yv(j)]';           
            xloc=RotMatrix'*(xglob-displ)+[dispx dispy]';
            if ( (xloc(IAXIS) > 0) && (xloc(IAXIS) < pi) && ...
                 (xloc(JAXIS) > 0) && (xloc(JAXIS) < pi) )
             v(i,j)= 0;
            else
             uloc(IAXIS)=-sin(t)*sin(k*xloc(IAXIS))^2*sin(2*k*xloc(JAXIS));
             uloc(JAXIS)= sin(t)*sin(2*k*xloc(IAXIS))*sin(k*xloc(JAXIS))^2;
             uglob= RotMatrix * uloc';          
             v(i,j)= uglob(JAXIS);
            end
        end
    end    
    
    % Cell centered velocities
    for j=1:ncell
        for i=1:ncell            
            xglob=[xc(i) yc(j)]';            
            xloc=RotMatrix'*(xglob-displ)+[dispx dispy]';
            if ( (xloc(IAXIS) > 0) && (xloc(IAXIS) < pi) && ...
                 (xloc(JAXIS) > 0) && (xloc(JAXIS) < pi) )
             uc(i,j)= 0;
             vc(i,j)= 0;
            else
             uloc(IAXIS)=-sin(t)*sin(k*xloc(IAXIS))^2*sin(2*k*xloc(JAXIS));
             uloc(JAXIS)= sin(t)*sin(2*k*xloc(IAXIS))*sin(k*xloc(JAXIS))^2;
             uglob= RotMatrix * uloc';
             uc(i,j)= uglob(IAXIS);
             vc(i,j)= uglob(JAXIS);
            end
        end
    end          
    
    % Pressure and Scalar:
    for j=1:ncell
        for i=1:ncell          
            xglob=[xc(i) yc(j)]';           
            xloc=RotMatrix'*(xglob-displ)+[dispx dispy]';
            
            P(i,j)= sin(t)/4*(2+cos(2*k*xloc(IAXIS)))*...
                             (2+cos(2*k*xloc(JAXIS)))-sin(t);
            
            Z(i,j)= A_Z/3*sin(t)*((1-cos(2*k*(xloc(IAXIS)-gam)))*...
                                  (1-cos(2*k*(xloc(JAXIS)-gam)))-1)+Z_mean;
            
            if ( (xloc(IAXIS) > 0) && (xloc(IAXIS) < pi) && ...
                 (xloc(JAXIS) > 0) && (xloc(JAXIS) < pi) )
               Z(i,j) = Z_mean;
            end
                        
            div(i,j) = (u(i+1,j)-u(i,j))/(xv(i+1)-xv(i)) + ...
                       (v(i,j+1)-v(i,j))/(yv(j+1)-yv(j));
        end
    end
    meanP = mean(mean(P))
    maxZ  = max(max(Z(1:ncell,1:ncell)));
    minZ  = min(min(Z(1:ncell,1:ncell))); 
    meanZ = mean(mean(Z));
    disp(['t=' num2str(t,'%6.3f') ...
          ', Min, Max, mean Z : ' num2str( minZ,'%8.5f') ...
                                ', ' num2str( maxZ,'%8.5f') ...
                                ', ' num2str(meanZ,'%13.10f')])
 
    % Scalar:
    contourf(xc(1:ncell)',yc(1:ncell)',Z(1:ncell,1:ncell)')
    % Rotated Cube:
    plot(xbd,ybd,'-w','LineWidth',5)
    % Velocities:
    quiver(xc(1:af:ncell),yc(1:af:ncell), ...
           uc(1:af:ncell,1:af:ncell)',vc(1:af:ncell,1:af:ncell)','w')
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
    title(['t=' num2str(t)],'FontSize',16)
    box on
    caxis([0.05 0.25])
    axis([xmin-0.1*Lx xmax+0.1*Lx ymin-0.1*Ly ymax+0.1*Ly]) 
    colorbar
    pause    
    if(it<nt)
      cla
    end     
end

meandiv=mean(mean((div)))
maxdiv=max(max(abs(div)))
