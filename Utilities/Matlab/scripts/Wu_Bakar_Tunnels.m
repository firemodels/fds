% McGrattan
% 1-22-2021
% Wu_Bakar_Tunnels.m
%
% This script replicates Fig. 4 of Wu and Bakar, Fire Safety Journal, 35 (2000) 363-390.

clear all
close all

outdir = '../../../out/Wu_Bakar_Tunnels/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Wu_Bakar_Tunnels/';

g = 9.81;
cp = 1.0;
Tinf = 293;
rho = 1.2;

label = {'A','B','C','D','E'};  % Tunnel names

Hbar = [0.176,0.250,0.333,0.400,0.238];       % Hydraulic tunnel diameters, 4*A/P
Q = [1.5,3.0,7.5,10.5,12.0,15.0,22.5,30.0];   % HRR (kW) of the propane burner

V = [0.43,0.46,0.48,0.48,0.48,0.48,0.48,0.48; ...
     0.39,0.48,0.56,0.59,0.60,0.60,0.60,0.60; ...
     0.37,0.45,0.54,0.57,0.59,0.60,0.62,0.65; ...
     0.34,0.40,0.50,0.54,0.56,0.59,0.65,0.65; ...
     0.44,0.54,0.60,0.60,0.60,0.60,0.60,0.60];  % Measured critical velocities (m/s)
L = zeros(5,8);

plot_style
Font_Interpreter='LaTeX';

for i=1:5

   % Find the back-layer length, L, defined as the distance from the upstream edge of the burner, 6.15 m, back to where
   % the back-layer temperature drops to 30 C, or 10 C above ambient.

   S = importdata([outdir,'Wu_Bakar_Tunnel_',label{i},'_line.csv'],',',2);

   clear indices
   for j=1:8
      indices = find( (S.data(:,j+1) > 30) & (S.data(:,j+1) < 50) );
      if size(indices)>0 L(i,j) = 6.15 - S.data(indices(1),1); end
   end

end

% Estimate the simulated critical velocity, Vmod, based on Li, Lei, and Ingason, FSJ 45 (2010) 361-370, Eq. 30.

for i=1:5
   for j=1:8
      Vmod(i,j) = V(i,j)*(1+L(i,j)/(18.5*Hbar(i)));
      Qstar(i,j) = Q(j)/(rho*cp*Tinf*sqrt(g)*Hbar(i)^2.5);
      Vstar(i,j) = Vmod(i,j)/sqrt(g*Hbar(i));
   end
end

% Difference between experimental critical velocity, V, and predicted, Vmod

Difference = 100*(Vmod-V)./V;

% Correlation of Wu and Bakar, Eq. 15-16, FSJ, 2000

Qcorr = [0.01:0.05:2];
for i=1:length(Qcorr)
   if Qcorr(i)< 0.2 Vcorr(i) = 0.4*0.2^(-1/3)*Qcorr(i)^(1/3); end
   if Qcorr(i)>=0.2 Vcorr(i) = 0.4; end
end

% Plot Qstar vs Vstar along with Qcorr vs Vcorr (Wu and Bakar, Fig. 4)

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
loglog(Qstar(1,:),Vstar(1,:),'kd'); hold on
loglog(Qstar(2,:),Vstar(2,:),'rs'); hold on
loglog(Qstar(3,:),Vstar(3,:),'m^'); hold on
loglog(Qstar(4,:),Vstar(4,:),'c*'); hold on
loglog(Qstar(5,:),Vstar(5,:),'go'); hold on
loglog(Qcorr,Vcorr,'k-'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
axis([0.001 10 0.1 1.0])
set(gca,'xticklabel',arrayfun(@(x) num2str(x),get(gca,'xtick'),'un',0))
set(gca,'yticklabel',arrayfun(@(y) num2str(y),get(gca,'ytick'),'un',0))
xlabel('Heat Release Rate, $Q^*$','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Critical Velocity, $V^*$','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
h = legend({'Tunnel A','Tunnel B','Tunnel C','Tunnel D','Tunnel E','Correlation'}, 'Location', 'NorthWest');
set(h,'Interpreter',Font_Interpreter)

% Add Git revision if file is available

Git_Filename = [outdir,'Wu_Bakar_Tunnel_A_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% Make the pdf figure

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Wu_Bakar_Critical_Velocity'])

display('Wu_Bakar_Tunnels completed successfully')
