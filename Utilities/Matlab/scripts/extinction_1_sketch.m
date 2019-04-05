% McGrattan
% 4-2-19
% extinction_1_sketch.m
%
% Makes figure for 'EXTINCTION 1' model for Tech Guide

close all
clear all

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

T1 = 0:600;
T2 = 600:1500;
X_O2_1 = 0.135*(1-(T1-20)/1427.);
X_O2_2 = 0.135*(1-(T2-20)/1427.);
XLegendStr{1}='constant volume';
XLegendStr{2}='quadratic';
XLegendStr{3}='user fan curve';
K(1)=plot(T1,X_O2_1,'k-','LineWidth',1); hold on
K(2)=plot(T2,X_O2_2,'k--','LineWidth',1);
K(3)=plot([600 600],[0 .08013],'k-','LineWidth',1);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 1500 0 0.2])
xlabel('Temperature (°C)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('Oxygen Volume Fraction','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
text(150,0.06,'No Burn','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
text(800,0.03,'Burn','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
text(700,0.15,'Burn','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
% print to pdf
plot_dir = '../../Manuals/FDS_Technical_Reference_Guide/FIGURES/';
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'extinction_1_sketch'])
