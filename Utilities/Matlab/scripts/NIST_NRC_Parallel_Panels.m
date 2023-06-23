% McGrattan
% 6-22-2023
% NIST_NRC_Parallel_Panels.m
%
% This script creates several different kinds of plots for NIST_NRC_Parallel_Panels cases.

clear all
close all

outdir = '../../../out/NIST_NRC_Parallel_Panels/';
expdir = '../../../exp/Submodules/macfp-db/Fire_Growth/NIST_Parallel_Panel/Experimental_Data/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/NIST_NRC_Parallel_Panels/';

plot_style

z = [10 20 30 50 75 100 140 180 220];
zm = linspace(2.5,243.5,50);
E = importdata([expdir,'PMMA_heatflux.csv'],',',2);
M = importdata([outdir,'PMMA_60_kW_1_cm_devc.csv'],',',2);

symbol = {'ko-' 'ro-' 'go-' 'co-' 'yo-' 'mo-' 'ro-'};

for i=[1 2 3 5 8 9 10]
   qdot{i} = E.data(i,2:10);
   err{i}  = E.data(i,11:19);
   errorbar(qdot{i},z,err{i},"horizontal","o")
   hold on
end

for i=[6 7 8 9 10 11 15]
   qdotm{i} = M.data(i,2:51);
   plot(qdotm{i},zm)
   hold on
end
   
xlabel('Heat Flux (kW/m²)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (cm)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
legend('120 kW','200 kW','300 kW','500 kW','1500 kW','2000 kW','2800 kW','Location','SouthEast','FontSize',10)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 150 0 250])

Git_Filename = [outdir,'PMMA_60_kW_1_cm_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'PMMA_Heat_Flux'])

hold off

