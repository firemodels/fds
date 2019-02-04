% McGrattan
% 10-21-2015
% scaling_tests.m
%
% Read _cpu.csv files for the MPI weak and strong scaling test cases

clear all
close all

plot_style

outdir = '../../../out/MPI_Scaling_Tests/FDS_Output_Files/';

M(1) = importdata([outdir,'strong_scaling_test_001_cpu.csv'],',',1);
M(2) = importdata([outdir,'strong_scaling_test_008_cpu.csv'],',',1);
M(3) = importdata([outdir,'strong_scaling_test_032_cpu.csv'],',',1);
M(4) = importdata([outdir,'strong_scaling_test_064_cpu.csv'],',',1);
M(5) = importdata([outdir,'strong_scaling_test_096_cpu.csv'],',',1);
M(6) = importdata([outdir,'strong_scaling_test_192_cpu.csv'],',',1);
M(7) = importdata([outdir,'strong_scaling_test_288_cpu.csv'],',',1);
M(8) = importdata([outdir,'strong_scaling_test_432_cpu.csv'],',',1);

r = [1 8 32 64 96 192 288 432];
r2 = [.1 8 32 64 96 192 432 1000];

for j=1:16
   for i=1:8
      t(i,j) = M(i).data(1,j)/M(1).data(1,16);
      t2(i) = 1./r2(i);
   end
end

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H1(1) = loglog(r2,2*t2,'k:'); hold on
H1(2) = loglog(r2,t2,'k:'); hold on
H1(3) = loglog(r2,t2/2,'k:'); hold on
H1(4) = loglog(r2,t2/4,'k:'); hold on
H1(5) = loglog(r2,t2/8,'k:'); hold on
H1(6) = loglog(r2,t2/16,'k:'); hold on
H1(7) = loglog(r2,t2/32,'k:'); hold on
H1(8) = loglog(r2,t2/64,'k:'); hold on
H1(9) = loglog(r2,t2/128,'k:'); hold on
H1(10) = loglog(r2,t2/256,'k:'); hold on

H(1) = loglog(r,t(:,16),'k-o');
H(2) = loglog(r,t(:,3),'r-o');
H(3) = loglog(r,t(:,4),'b-o');
H(4) = loglog(r,t(:,5),'m-o');
H(5) = loglog(r,t(:,6),'c-o');
H(6) = loglog(r,t(:,12),'g-o');
H(7) = loglog(r,t(:,10),'y-o');
H(8) = loglog(r,t(:,2),'k-s');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('MPI Processes','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Relative Wall Clock Time','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Min_Ind = 1.0;
Max_Ind = 1000;
Min_Dep = 0.00001;
Max_Dep = 1.;
axis([Min_Ind Max_Ind Min_Dep Max_Dep])
set(gca,'XTickLabel',num2str(get(gca,'XTick')'))
yticks([0.000001 0.00001 0.0001 0.001 0.01 0.1 1]);
Title_Position(1) = 0.40;
Title_Position(2) = 0.95;
X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));
text(X_Title_Position,Y_Title_Position,'Strong Scaling Test','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
legend_handle = legend(H,'Total','DIVG','MASS','VELO','PRES','COMM','RADI','MAIN','Location','NorthEast');
set(legend_handle,'Interpreter',Font_Interpreter);
set(legend_handle,'Fontsize',8);

git_file = [outdir,'strong_scaling_test_288_git.txt'];
addverstr(gca,git_file,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/strong_scaling_test'])

close all
clear all
plot_style

outdir = '../../../out/MPI_Scaling_Tests/FDS_Output_Files/';

M(1)  = importdata([outdir,'weak_scaling_test_001_cpu.csv'],',',1);
M(2)  = importdata([outdir,'weak_scaling_test_002_cpu.csv'],',',1);
M(3)  = importdata([outdir,'weak_scaling_test_004_cpu.csv'],',',1);
M(4)  = importdata([outdir,'weak_scaling_test_008_cpu.csv'],',',1);
M(5)  = importdata([outdir,'weak_scaling_test_016_cpu.csv'],',',1);
M(6)  = importdata([outdir,'weak_scaling_test_032_cpu.csv'],',',1);
M(7)  = importdata([outdir,'weak_scaling_test_064_cpu.csv'],',',1);
M(8)  = importdata([outdir,'weak_scaling_test_128_cpu.csv'],',',1);
M(9)  = importdata([outdir,'weak_scaling_test_192_cpu.csv'],',',1);
M(10) = importdata([outdir,'weak_scaling_test_288_cpu.csv'],',',1);
M(11) = importdata([outdir,'weak_scaling_test_432_cpu.csv'],',',1);

r = [1 2 4 8 16 32 64 128 192 288 432];

for i=1:11
   t(i) = M(1).data(1,16)/M(i).data(1,16);
   t2(i) = 1.;
end

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1) = semilogx(r,t,'ko'); hold on
H(2) = semilogx(r,t2,'k--');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('MPI Processes','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Efficiency','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Min_Ind = 1.0;
Max_Ind = 1000;
Min_Dep = 0.0;
Max_Dep = 1.2;
axis([Min_Ind Max_Ind Min_Dep Max_Dep])
set(gca,'XTickLabel',num2str(get(gca,'XTick')'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')'))
Title_Position(1) = 0.60;
Title_Position(2) = 0.90;
X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
text(X_Title_Position,Y_Title_Position,'Weak Scaling Test','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
legend(H,'FDS','Ideal','Location','Southwest')

git_file = [outdir,'weak_scaling_test_288_git.txt'];
addverstr(gca,git_file,'semilogx')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/weak_scaling_test'])

