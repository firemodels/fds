% McDermott and Floyd
% 10-3-12
% fan_curve.m
%
% Makes figure for HVAC Fan Parameters section of the FDS Users Guide

close all
clear all

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

vdot_max = 10;
dp_max = 500;

dp = -1000:1000;

vdot = vdot_max*sign(dp_max-dp).*sqrt(abs(dp-dp_max)/dp_max);
vdot1 = 10;
XLegendStr{1}='constant volume';
XLegendStr{2}='quadratic';
XLegendStr{3}='user fan curve';
K(1)=plot(vdot1*ones(1,length(dp)),dp,'r-','LineWidth',2); hold on
K(2)=plot(vdot,dp,'k-','LineWidth',2);
i=0;
for dp=-1000:200:1000
    i=i+1;
    rampx(i) = vdot_max*sign(dp_max-dp).*sqrt(abs(dp-dp_max)/dp_max);
    rampy(i) = dp;
end
K(3)=plot(rampx,rampy,'b-','LineWidth',2);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'XGrid','on')
set(gca,'YGrid','on')
axis([-10 20 -1000 1000])
xlabel('Volume Flow Rate (m^3/s)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('Static Pressure (Pa)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
lh=legend(K,XLegendStr,'Location','Southwest');
set(lh,'FontSize',Key_Font_Size)
% print to pdf
plot_dir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'fan_curve'])
