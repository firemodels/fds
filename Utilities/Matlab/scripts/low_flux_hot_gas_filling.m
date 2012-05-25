% McDermott
% 9-24-10
% low_flux_hot_gas_filling.m

close all
clear all

addpath('../../../Verification/Flowfields')

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

filename = 'low_flux_hot_gas_filling_mass.csv';
H = plot_mass(filename,.0001,6,2,4);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Accumulated Mass (kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 60 0 .036])
text(5,.032,'Low Flux Hot Gas Filling','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
hh=legend(H,'Specified','FDS','Location','Southeast');
legend boxon
set(hh,'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size)

% add SVN if file is available

SVN_Filename = 'low_flux_hot_gas_filling_svn.txt';
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

plotdir = '../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
Plot_Filename = 'low_flux_hot_gas_filling';

set(gcf,'Visible','on');
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
display('Printing plot low_flux_hot_gas_filling.pdf...')
print(gcf,'-dpdf',[plotdir,Plot_Filename])