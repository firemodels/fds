% mesh_quad_stretch : One dimensional grid quadratic stretching.
%
% This script provides a set of FDS &TRNX, &TRNY or &TRNZ grid stretching 
% input lines. It is of use on case where the cell size is required to 
% increase with position, and the cell size DXS on the low boundary needs
% to be prescribed.
%
% The known paramaters to be prescribed are:
%   - the length of mesh block Lx, and start of block xs.
%   - the uniform grid size DXS of the first m cells,
%   - and the total number of cells in the mesh Nx.
%
% We start by defining:
% xf = xs + sum_{i=1}^{Nx} DXi
% where the boundary positions of the grid are xs, xf, Lx=xf-xs.
%
% => Lx = m*DXS + sum_{i=m+1}^{Nx} DXi
%
% Here we assume DXi = a DXi-1, where the parameter a is a constant.
%
% In the previous:
% => Lx = m*DXS + DXS*sum_{i=m+1}^{Nx} a^(i-m)
%
% Therefore, the polynomial:
% (m-Lx/DXS) + sum_{i=m+1}^{Nx} a^(i-m) = 0 is solved for a, around a0=1.
%
% In order to get increasing cell sizes: DXS < Lx/Nx.
%
% -------------------------------------------------------------------------
close all
clear all
clc

%% Directory where to write TRN.dat file: 
basedir='MY DIRECTORY';

% Write out TRNX, TRNY, TRNZ?
TRN='TRNZ';
TRN_ID='''MY TRANSFORM''';

% Mesh size:
xs = 1; % Mesh block is defined from xs=1 to xf=xs+Lx=2.
Lx = 1;

% Low side number m of cells with small uniform size DXS:
m  =     4;
DXS= 0.025;

% Total number of cells:
Nx =   20;

%% Solve for a:
disp('Solving for stretching parameter..')

% Define tolerance:
etol_rel = 1.e-10;
max_iter = 1000;

% Newton iteration:
t=cputime;
an = 1.;
for n=1:max_iter
    % Compute f, fp:
    suma = 0;;
    sumap= 0.;
    for i=m+1:Nx
        suma = suma  + an^(i-m);
        sumap= sumap + (i-m)*an^(i-m-1);
    end
    f_an = (m-Lx/DXS) + suma;
    fp_an= sumap;
    
    % New value of a:
    an1 = an - f_an/fp_an;
    
    % Convergence test:
    if(abs(an1-an)/an < etol_rel)
        disp(['Iter ' num2str(n,'%4.4d') ...
             ', Convergence found, stretch factor a=' num2str(an1) ... 
             ', relative error=' num2str(abs(an1-an)/an) '.'])
        break
    end
    if(mod(n,ceil(max_iter/1000)) == 0)
        disp(['Iter ' num2str(n,'%4.4d') ...
             ', Relative error=' num2str(abs(an1-an)/an,'%18.12f') '.'])
    end
    % Update an:
    an=an1;
end
disp(['Time taken :' num2str(cputime-t) ' sec.'])

%% Make figure:
figure
subplot(1,2,1)
hold on
plot([0.1*Lx 0.1*Lx],[xs xs+Lx],'k')
plot([0.3*Lx 0.3*Lx],[xs xs+Lx],'k')
plot([0.1*Lx 0.3*Lx],[xs xs],'k')
xpl(1) = xs; xipl(1) = xs;
for i=1:m
   plot([0.1*Lx 0.3*Lx],[xs+i*DXS xs+i*DXS],'k')
   xpl(i+1) =xs+i*DXS;
   xipl(i+1)=xs+i*Lx/Nx;
end
for i=m+1:Nx
    dx = DXS;
    sdx= 0.;
    for k=m+1:i
        dx = an1*dx;
        sdx= sdx + dx;
    end
    x = m*DXS + sdx;
    plot([0.1*Lx 0.3*Lx],[xs+x xs+x],'k')
    xpl(i+1) =xs+x;
    xipl(i+1)=xs+i*Lx/Nx;
end
xpl(Nx+1) =xs+Lx;
xipl(Nx+1)=xs+Lx;
box on
ylabel('Stretching direction','FontSize',16)
set(gca,'XTick',[],'FontSize',14)
axis([0 0.4*Lx xs-0.1*Lx xs+1.1*Lx]);

subplot(1,2,2)
plot(xipl,xpl,'+k','LineWidth',2)
xlabel('\xi','FontSize',16)
ylabel('Physical coordinate','FontSize',16)
set(gca,'FontSize',14)
axis equal
box on
grid on

%% Write transformation input lines:
disp(' ')
fprintf(['Writing stretching input file...\n'])
fprintf([basedir 'TRN.dat\n'])

[fid]=fopen([basedir 'TRN.dat'],'w');
for i=1:m
   % First m cells:
   fprintf(fid,['&' TRN ' ID=' TRN_ID ', CC=%20.15f, PC=%20.15f /\n'],...
            xs+i*Lx/Nx,xs+i*DXS);
end
for i=m+1:Nx-1
    dx = DXS;
    sdx= 0.;
    for k=m+1:i
        dx = an1*dx;
        sdx= sdx + dx;
    end
    x = m*DXS + sdx;
    % Stretched cells:
    fprintf(fid,['&' TRN ' ID=' TRN_ID ', CC=%20.15f, PC=%20.15f /\n'],...
            xs+i*Lx/Nx,xs+x);
end
fclose(fid);
fprintf(['Done.\n'])