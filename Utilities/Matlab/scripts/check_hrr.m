% McDermott
% 9-20-10
% check_hrr.m

close all
clear all

% set the plot style parameters

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

DDIR = '../../Validation/Flame_Height/FDS_Output_Files/';

M = csvread([DDIR,'box_height.csv'],1,0);
Qs = M(:,1);
Q = M(:,2);
RI = {'_RI=05','','_RI=20'};
QI = {'p1','p2','p5','1','2','5','10','20','50','100','200','500','1000','2000','5000','10000'};

K(4)=loglog(Qs,Q,'k-'); hold on

for j=1:length(Qs)
    for i=1:3
        filename = ['Qs=',QI{j},RI{i},'_hrr.csv'];
        A = csvread([DDIR,filename],2,0);
        n = length(A(:,1));
        Q_fds = mean(A(round(n/2):n,2));
        if (i==1); K(1)=loglog(Qs(j),Q_fds,'ksq'); end
        if (i==2); K(2)=loglog(Qs(j),Q_fds,'r^'); end
        if (i==3); K(3)=loglog(Qs(j),Q_fds,'go'); end
    end
end

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Q*','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('FDS Q (kW)','Interpreter','LaTeX','FontSize',Label_Font_Size)
text(.3,2e7,'Flame Height Heat Release Verification','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter','LaTeX')
hh=legend(K,'$D^*/\delta x=5$','$D^*/\delta x=10$','$D^*/\delta x=20$','correct','Location','Southeast');
set(hh,'Interpreter','LaTeX','FontSize',Key_Font_Size)

% print to pdf

plotdir = '../../Manuals/FDS_5_Validation_Guide/FIGURES/Heskestad/';
Plot_Filename = 'Flame_Height_check_hrr';

set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
display(['Printing plot Flame_Height_check_hrr.pdf ...'])
print(gcf,'-dpdf',[plotdir,Plot_Filename])



