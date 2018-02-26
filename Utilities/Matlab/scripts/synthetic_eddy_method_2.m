% McGrattan
% 2-26-2018
% synthetic_eddy_method_2.m
% This script should probably be appended to synthetic_eddy_method.m

close all
clear all

plot_style
Font_Interpreter = 'LaTeX';

datadir='../../Verification/Turbulence/';
plotdir='../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

symbol = {'b*','rs','bx','go','r^'};
leg    = {'DeGraff and Eaton (2000)','Lund et al. (1998)','Rai and Moin (1993)','Spalart (1988)','Wu et al. (2009)'};

% --------
% z+ vs u+
% --------
figure
set(gca,'Units',Plot_Units)
set(gca,'Fontname',Font_Name)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([datadir,'BL_Stats_UplusZplus.csv'],',',1);
label = {'u+_Eaton','u+_Lund','u+_Rai','u+_Spalart','u+_Wu'};

for i=1:5
   k = find(strcmp(M.colheaders,label(i)));
   z = M.data(:,k-1);
   u = M.data(:,k);
   H(i)=semilogx(z,u,symbol{i}); hold on
 % clear k,z,u
end

axis([1 1e7 0 50])
set(gca,'Fontname',Font_Name)
ylabel('$u^+$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z^+$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

lh = legend(H,leg{1},leg{2},leg{3},leg{4},leg{5},'Location','SouthEast');
set(lh,'FontSize',Key_Font_Size,'Fontname',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[plotdir,'BL_Stats_UplusZplus'])
clear H lh

% ---------------
% z/delta vs u'u'
% ---------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([datadir,'BL_Stats_UprimeUprime.csv'],',',1);
label = {'up_up_Eaton','up_up_Lund','up_up_Rai','up_up_Spalart','up_up_Wu'};

for i=1:5
   k = find(strcmp(M.colheaders,label(i)));
   z = M.data(:,k-1);
   u = M.data(:,k);
   H(i)=plot(z,u,symbol{i}); hold on
 % clear k,z,u
end

axis([0 1.25 0 15])
set(gca,'Fontname',Font_Name)
ylabel('$\overline{u^\prime u^\prime}/u_\tau^2$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z/\delta$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

lh = legend(H,leg{1},leg{2},leg{3},leg{4},leg{5},'Location','NorthEast');
set(lh,'FontSize',Key_Font_Size,'Fontname',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[plotdir,'BL_Stats_UprimeUprime'])
clear H lh M label symbol

% ---------------
% z/delta vs u'w'
% ---------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([datadir,'BL_Stats_UprimeWprime.csv'],',',1);
label = {'up_wp_Lund','up_wp_Spalart','up_wp_Wu'};
symbol = {'rs','go','r^'};

for i=1:3
   k = find(strcmp(M.colheaders,label(i)));
   z = M.data(:,k-1);
   u = M.data(:,k);
   H(i)=plot(z,u,symbol{i}); hold on
 % clear k,z,u
end

axis([0 1.25 -2 2])
set(gca,'Fontname',Font_Name)
ylabel('$\overline{u^\prime w^\prime}/u_\tau^2$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z/\delta$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

lh = legend(H,leg{2},leg{4},leg{5},'Location','NorthEast');
set(lh,'FontSize',Key_Font_Size,'Fontname',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[plotdir,'BL_Stats_UprimeWprime'])
clear H lh M label symbol


% ---------------
% z/delta vs v'v'
% ---------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([datadir,'BL_Stats_VprimeVprime.csv'],',',1);
label = {'vp_vp_Eaton','vp_vp_Lund','vp_vp_Spalart','vp_vp_Wu'};
symbol = {'b*','rs','go','r^'};

for i=1:4
   k = find(strcmp(M.colheaders,label(i)));
   z = M.data(:,k-1);
   u = M.data(:,k);
   H(i)=plot(z,u,symbol{i}); hold on
end

axis([0 1.25 0 6])
set(gca,'Fontname',Font_Name)
ylabel('$\overline{v^\prime v^\prime}/u_\tau^2$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z/\delta$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

lh = legend(H,leg{1},leg{2},leg{4},leg{5},'Location','NorthEast');
set(lh,'FontSize',Key_Font_Size,'Fontname',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[plotdir,'BL_Stats_VprimeVprime'])
clear H lh M label symbol

% ---------------
% z/delta vs w'w'
% ---------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([datadir,'BL_Stats_WprimeWprime.csv'],',',1);
label = {'wp_wp_Eaton','wp_wp_Lund','wp_wp_Spalart','wp_wp_Wu'};
symbol = {'b*','rs','go','r^'};

for i=1:4
   k = find(strcmp(M.colheaders,label(i)));
   z = M.data(:,k-1);
   u = M.data(:,k);
   H(i)=plot(z,u,symbol{i}); hold on
end

axis([0 1.25 0 4])
set(gca,'Fontname',Font_Name)
ylabel('$\overline{w^\prime w^\prime}/u_\tau^2$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z/\delta$','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

lh = legend(H,leg{1},leg{2},leg{4},leg{5},'Location','NorthEast');
set(lh,'FontSize',Key_Font_Size,'Fontname',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[plotdir,'BL_Stats_WprimeWprime'])
clear H lh M label symbol

