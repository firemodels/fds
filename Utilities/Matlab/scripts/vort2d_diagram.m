%-------------------------------------------------------------------------%
% FILE: vort2d_diagram.m
% AUTHOR: Max Gould
% LAST DATE EDITED: 06 August 2012
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Generate flow diagram for the 2D Vortex section
% of the FDS Verification Guide
%
% Dependencies:
%    ../scripts/arrows.m
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Define Parameters
%-------------------------------------------------------------------------%

% Define basic flow parameters
U0 = 0.7;
Rc = 0.01556;
BigGamma = 0.0359157;

% Define plot parameters
Min  = -0.05;
Max  =  0.05;
Step =  0.005;

% Define plot region
[X,Z] = meshgrid(Min:Step:Max);

%-------------------------------------------------------------------------%
% Setup Flow Field
%-------------------------------------------------------------------------%

% Define flow field
Psi = BigGamma.*exp(-(X.^2+Z.^2)./(2.*Rc.^2));

% Calculate gradient and define velocity components
[DX,DZ] = gradient(Psi,0.005,0.005);
U = U0+DZ;
W = -DX;

%-------------------------------------------------------------------------%
% Plot Vector Field
%-------------------------------------------------------------------------%

% Define standard plot-style parameters
plot_style

% Plot vector field
quiver(X,Z,U,W,2,'Color',[0.3 0.7 1],'ShowArrowHead','off');
hold on

% Hold plot for arrow heads
arrows(X,Z,U,W,0.32,0.25);
hold off

%-------------------------------------------------------------------------%
% Plot Arrow Heads Without arrows.m
%-------------------------------------------------------------------------%
%     %-----------------------------------------------------------------%
%     % Plot Arrow Heads
%     %-----------------------------------------------------------------%
%
%     % Define arrow-head dimension and scaling parameters
%     Alpha = 0.32;      %primary scaling factor (arrow size)
%     Beta  = 0.25;      %secondary scaling factor (width to length ratio)
%     Scale = 1.6*Step;  %scales arrow-head placement with vector lines
% 
%     % Convert arrays to vectors and scale
%     X=X(:).';
%     Z=Z(:).';
%     U=Scale*U(:).';
%     W=Scale*W(:).';
% 
%     % Allocate space for arrays
%     hu = zeros(length(U),4);
%     hw = zeros(length(U),4);
% 
%     % Define arrow geometry and plot arrows
%     for i = 1:length(U)
%         hu(i,1) = X(i)+U(i)-Alpha*(U(i)+Beta*(W(i)+eps));
%         hu(i,2) = X(i)+U(i);
%         hu(i,3) = X(i)+U(i)-Alpha*(U(i)-Beta*(W(i)+eps));
%         hu(i,4) = X(i)+U(i)-Alpha*(U(i)+Beta*(W(i)+eps));
%         hw(i,1) = Z(i)+W(i)-Alpha*(W(i)-Beta*(U(i)+eps));
%         hw(i,2) = Z(i)+W(i);
%         hw(i,3) = Z(i)+W(i)-Alpha*(W(i)+Beta*(U(i)+eps));
%         hw(i,4) = Z(i)+W(i)-Alpha*(W(i)-Beta*(U(i)+eps));
%         h2 = fill(hu(i,:),hw(i,:),[0.3 0.7 1]);
%     end
% 
% % Plot arrow heads with vector lines
% hold off
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Produce Plot
%-------------------------------------------------------------------------%

% Define output location
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Define plot range
axis([-0.05,0.05,-0.05,0.05]);

set(gca,'XTick',-0.04:0.02:0.04);
set(gca,'YTick',-0.04:0.02:0.04);
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Width])
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Width]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Width]);

% Save plot to file
print(gcf,'-dpdf',[plot_dir,'vort2d_diagram']);

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%