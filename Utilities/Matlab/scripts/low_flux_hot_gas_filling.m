% McDermott
% 9-24-10
% low_flux_hot_gas_filling.m

close all
clear all

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

filename = '../../../Verification/Flowfields/low_flux_hot_gas_filling_mass.csv';
H = plot_mass(filename,.0001,6,2,3);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Time (s)','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('Accumulated Mass (kg)','Interpreter','LaTeX','FontSize',Label_Font_Size)
axis([0 60 0 .036])
text(5,.032,'Low Flux Hot Gas Filling','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter','LaTeX')
hh=legend(H,'Specified','FDS','Location','Southeast');
legend boxon
set(hh,'Interpreter','LaTeX','FontSize',Key_Font_Size)

% print to pdf

plotdir = '../../../Manuals/FDS_5_Verification_Guide/SCRIPT_FIGURES/';
Plot_Filename = 'low_flux_hot_gas_filling';

set(gcf,'Visible','on');
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
display('Printing plot low_flux_hot_gas_filling.pdf...')
print(gcf,'-dpdf',[plotdir,Plot_Filename])