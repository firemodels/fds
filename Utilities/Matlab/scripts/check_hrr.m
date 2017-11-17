% McDermott
% 9-20-10
% check_hrr.m

close all
clear all

% set the plot style parameters

plot_style

outdir = '../../../out/Heskestad_Flame_Height/FDS_Output_Files/';
M = csvread([outdir,'box_height.csv'],1,0);
Qs = M(:,1);
Q = M(:,2);
RI = {'_RI=05','_RI=10','_RI=20'};
QI = {'p1','p2','p5','1','2','5','10','20','50','100','200','500','1000','2000','5000','10000'};

K(4)=loglog(Qs,Q,'k-'); hold on

for j=1:length(Qs)
    for i=1:3
        filename = ['Qs=',QI{j},RI{i},'_hrr.csv'];
        A = csvread([outdir,filename],2,0);
        n = length(A(:,1));
        Q_fds = mean(A(round(n/2):n,2));
        if (i==1); K(1)=loglog(Qs(j),Q_fds,'ksq'); end
        if (i==2); K(2)=loglog(Qs(j),Q_fds,'r^'); end
        if (i==3); K(3)=loglog(Qs(j),Q_fds,'go'); end
    end
end

axis([.05 10^4 10^2 10^8])

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

xlabel('{\it Q*}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Heat Release Rate (kW)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(.1,2e7,'Flame Height Heat Release Verification','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
hh=legend(K,'{\itD}*/{\it\deltax}=5','{\itD}*/{\it\deltax}=10','{\itD}*/{\it\deltax}=20','correct','Location','Southeast');
set(hh,'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size)

% add version if file is available

Git_Filename = [outdir,'Qs=1_RI=05_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf

plotdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Heskestad/';
Plot_Filename = 'Flame_Height_check_hrr';

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
display(['Printing plot Flame_Height_check_hrr.pdf ...'])
print(gcf,'-dpdf',[plotdir,Plot_Filename])



