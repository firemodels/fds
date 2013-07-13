% Hostikka
% 8-11-2009
% birch_tga.m

close all
clear all

addpath('../../Verification/Pyrolysis')

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)

TGA_20_N2 = csvread('birch_tga_20_exp.csv',34,0);
TGA_2_N2 = csvread('birch_tga_2_exp.csv',34,0);

skip_case = 0;
if ~exist('birch_tga_1step_2_devc.csv')
    display('Error: File birch_tga_1step_2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('birch_tga_1step_20_devc.csv')
    display('Error: File birch_tga_1step_20_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if skip_case
    return
end

M_2 = csvread('birch_tga_1step_2_devc.csv',2,0);
M_20 = csvread('birch_tga_1step_20_devc.csv',2,0);

h=plot(TGA_2_N2(:,1),TGA_2_N2(:,2)/100,TGA_20_N2(:,1),TGA_20_N2(:,2)/100,M_2(:,3),M_2(:,2)/M_2(1,2),M_20(:,3),M_20(:,2)/M_20(1,2));
hold on
set(gca,'FontName',Font_Name)
axis([0 800 0 1.1])
xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(h([1 3]),'Color','b')
set(h([2 4]),'Color','r')
set(h([3 4]),'LineStyle','--')
legend('Exp. 2 K/min','Exp. 20 K/min','FDS 2 K/min','FDS 20 K/min')

% add SVN if file is available

SVN_Filename = ['birch_tga_1step_2_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/birch_tga')
