% McGrattan
% 10-21-2015
% scaling_tests.m
%
% Read _cpu.csv files for the MPI weak and strong scaling test cases

clear all
close all

FDS_Output_Files = '../../Validation/MPI_Scaling_Tests/FDS_Output_Files/';

M(1) = importdata([FDS_Output_Files,'strong_scaling_test_001_cpu.csv'],',',1);
M(2) = importdata([FDS_Output_Files,'strong_scaling_test_008_cpu.csv'],',',1);
M(3) = importdata([FDS_Output_Files,'strong_scaling_test_032_cpu.csv'],',',1);
M(4) = importdata([FDS_Output_Files,'strong_scaling_test_064_cpu.csv'],',',1);
M(5) = importdata([FDS_Output_Files,'strong_scaling_test_096_cpu.csv'],',',1);
M(6) = importdata([FDS_Output_Files,'strong_scaling_test_192_cpu.csv'],',',1);
M(7) = importdata([FDS_Output_Files,'strong_scaling_test_288_cpu.csv'],',',1);

r = [1 8 32 64 96 192 288];

for i=1:7
   t(i) = M(i).data(1,15)/M(1).data(1,15);
   t2(i) = 1./r(i);
end

plot_style

H(1) = loglog(r,t,'ko'); hold on
H(2) = loglog(r,t2,'k--'); 
 
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('MPI Processes','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Relative CPU Time','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Min_Ind = 1;
Max_Ind = 300;
Min_Dep = 0.001;
Max_Dep = 1.;
axis([Min_Ind Max_Ind Min_Dep Max_Dep])
Title_Position(1) = 0.60;
Title_Position(2) = 0.90;
X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));
text(X_Title_Position,Y_Title_Position,'Strong Scaling Test','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
legend(H,'FDS','Ideal','Location','Southwest')

git_file = [FDS_Output_Files,'strong_scaling_test_288_git.txt'];
addverstr(gca,git_file,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/strong_scaling_test'])

clear all
close all

FDS_Output_Files = '../../Validation/MPI_Scaling_Tests/FDS_Output_Files/';

M(1)  = importdata([FDS_Output_Files,'weak_scaling_test_001_cpu.csv'],',',1);
M(2)  = importdata([FDS_Output_Files,'weak_scaling_test_002_cpu.csv'],',',1);
M(3)  = importdata([FDS_Output_Files,'weak_scaling_test_004_cpu.csv'],',',1);
M(4)  = importdata([FDS_Output_Files,'weak_scaling_test_008_cpu.csv'],',',1);
M(5)  = importdata([FDS_Output_Files,'weak_scaling_test_016_cpu.csv'],',',1);
M(6)  = importdata([FDS_Output_Files,'weak_scaling_test_032_cpu.csv'],',',1);
M(7)  = importdata([FDS_Output_Files,'weak_scaling_test_064_cpu.csv'],',',1);
M(8)  = importdata([FDS_Output_Files,'weak_scaling_test_128_cpu.csv'],',',1);
M(9)  = importdata([FDS_Output_Files,'weak_scaling_test_192_cpu.csv'],',',1);
M(10) = importdata([FDS_Output_Files,'weak_scaling_test_288_cpu.csv'],',',1);

r = [1 2 4 8 16 32 64 128 192 288];

for i=1:10
   t(i) = M(1).data(1,15)/M(i).data(1,15);
   t2(i) = 1.;
end

plot_style

H(1) = semilogx(r,t,'ko'); hold on
H(2) = semilogx(r,t2,'k--');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('MPI Processes','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Efficiency','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Min_Ind = 1;
Max_Ind = 300;
Min_Dep = 0.0;
Max_Dep = 1.2;
axis([Min_Ind Max_Ind Min_Dep Max_Dep])
Title_Position(1) = 0.60;
Title_Position(2) = 0.90;
X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
text(X_Title_Position,Y_Title_Position,'Weak Scaling Test','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
legend(H,'FDS','Ideal','Location','Southwest')

git_file = [FDS_Output_Files,'weak_scaling_test_288_git.txt'];
addverstr(gca,git_file,'semilogx')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/weak_scaling_test'])

