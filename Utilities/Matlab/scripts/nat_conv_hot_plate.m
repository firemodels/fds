% McGrattan
% September 1, 2023
% nat_conv_hot_plate.m

close all
clear all

plot_style

results_dir = ['../../../out/Convection/'];

% calculations below were used for input file setup

g  = 9.80665;
S  = [1.0  0.1  5.0  10 .05];
T1 = [503  503  503 323 323];
T2 = 293;
MW = 28.85476; % FDS 'LJ AIR'
P0 = 101325;
mu = 1.8216e-5;
cp = 1000;
k  = 0.018216; % for Pr=1 fluid

RAYLEIGH_DOWN = logspace(4, 9,100);
RAYLEIGH_UP   = logspace(4,11,100);
for i=1:length(RAYLEIGH_DOWN)
   NUSSELT_DOWN(i)=0.52*RAYLEIGH_DOWN(i)^0.2;  % Incropera and DeWitt, 7th, Eq. 9.32
end
for i=1:length(RAYLEIGH_UP)
   if RAYLEIGH_UP(i)<1e7; 
      NUSSELT_UP(i)=0.54*RAYLEIGH_UP(i)^0.25;  % Incropera and DeWitt, 7th, Eq. 9.30
   else
      NUSSELT_UP(i)=0.15*RAYLEIGH_UP(i)^0.333; % Incropera and DeWitt, 7th, Eq. 9.31
   end
end

figure
set(gcf,'Visible',Figure_Visibility);
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

marker_handle(1)=loglog(RAYLEIGH_DOWN,NUSSELT_DOWN,'k-'); hold on
marker_handle(2)=loglog(RAYLEIGH_UP  ,NUSSELT_UP  ,'k--'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
axis([1e0 1e15 1 1000])
ytick = [1 2 5 10 20 50 100 200 500 1000];
yticks(ytick);
xlabel('Rayleigh Number','FontSize',Label_Font_Size)
ylabel('Nusselt Number','FontSize',Label_Font_Size)

% FDS results

casename={...
'nat_conv_hot_plate_1',...
'nat_conv_hot_plate_2',...
'nat_conv_hot_plate_3',...
'nat_conv_hot_plate_4',...
'nat_conv_hot_plate_5',...
};

marker_style = {'gsq','rsq','ksq','go','ro','ko'};
res = {'8','16','32'};

for j=1:length(res)
    for i=1:length(S)

        M = importdata([results_dir,casename{i},'_',res{j},'_devc.csv']);

        % check for steady state
        Qdot_down = 1000*M.data(end,2);
        Qdot_up = 1000*M.data(end,3);
        A = S(i)^2;
        Tm = 0.5*(T1(i)+T2);
        beta = 1/Tm;
        rho = P0*MW/(8341.5*Tm);
        alpha = k/(rho*cp);
        nu = mu/rho;
        L = S(i)/4;
        Ra_FDS = (g*beta*(T1(i)-T2)*L^3)/(alpha*nu);
        Nu_FDS_down = (Qdot_down/A)*(L/k)/(T1(i)-T2);
        Nu_FDS_up   = (Qdot_up  /A)*(L/k)/(T1(i)-T2);

        marker_handle(j+2)=plot(Ra_FDS,Nu_FDS_down,marker_style{j},'MarkerSize',8);
        marker_handle(j+5)=plot(Ra_FDS,Nu_FDS_up,marker_style{j+3},'MarkerSize',8);
    end
end

lh=legend(marker_handle,'Correlation; Down','Correlation; Up','Down ($S/\delta x=8$)','Down ($S/\delta x=16$)','Down ($S/\delta x=32$)',...
  'Up ($S/\delta x=8$)','Up ($S/\delta x=16$)','Up ($S/\delta x=32$)','Location','Northwest','Interpreter','LaTeX');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'nat_conv_hot_plate_1_8_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/nat_conv_hot_plate');


