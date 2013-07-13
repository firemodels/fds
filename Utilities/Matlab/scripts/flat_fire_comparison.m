% Trettel
% flat fire comparison
% note that fire here is a verb, as in firing a gun, as this approximate solution is for ballistic trajectories

close all
clear all

V_0 = 400;
h = 8;
g = 9.8;
C_d = 0.2;
rho_a = 1.2;
rho_d = 1000;
D = 5e-3;

K = 3 * rho_a * C_d / (4 * rho_d * D);

dt = 0.05;
tend = 1.65;

tvec = 0:dt:tend;

xexact = log(V_0 * K * tvec + 1) / K;
yexact = h...
    + (g / (2 * (K * V_0)^2)) * log(V_0 * K * tvec + 1)...
    - g * tvec .^ 2 / 4 ...
    - g * tvec ./ (2 * K * V_0);
uexact = V_0 ./ (V_0 * K * tvec + 1);
vexact = g ./ (2 * K * V_0 * (K * V_0 * tvec + 1))...
    - g * tvec ./ 2 ...
    - g / (2 * K * V_0);

repository = '../../Verification/Sprinklers_and_Sprays/';

if ~exist([repository, 'flat_fire.prt5'])
    display(['Error: File ' [repository, 'flat_fire.prt5'] ' does not exist. Skipping case.'])
    return
end

[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'flat_fire.prt5'],'real*4');

H(1) = plot(xexact, yexact, '-k');
hold on
H(2) = plot(XP, ZP, 'ok');
hold off
plot_style
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)
set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)
xlabel('{\it x}-position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it z}-position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
h = legend(H, 'Flat fire', 'FDS', 'Location', 'Southwest');
set(h,'Interpreter', Font_Interpreter)
set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);
SVN_Filename = [repository, 'flat_fire_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end
display('Printing plot flat_fire_trajectory.pdf...')
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/flat_fire_trajectory')

close all
H(1) = plot(tvec, xexact, '-k');
hold on
H(2) = plot(STIME, XP, 'ok');
hold off
plot_style
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)
set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)
xlabel('Time (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it x}-position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
h = legend(H, 'Flat fire', 'FDS', 'Location', 'Southeast');
set(h,'Interpreter', Font_Interpreter)
set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);
SVN_Filename = [repository, 'flat_fire_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end
display('Printing plot flat_fire_x.pdf...')
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/flat_fire_x')

close all
H(1) = plot(tvec, yexact, '-k');
hold on
H(2) = plot(STIME, ZP, 'ok');
hold off
plot_style
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)
set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)
xlabel('Time (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it z}-position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
h = legend(H, 'Flat fire', 'FDS', 'Location', 'Southwest');
set(h,'Interpreter', Font_Interpreter)
set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);
SVN_Filename = [repository, 'flat_fire_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end
display('Printing plot flat_fire_z.pdf...')
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/flat_fire_z')

close all
H(1) = plot(tvec, uexact, '-k');
hold on
H(2) = plot(STIME, QP(:,:,1,1), 'ok');
hold off
plot_style
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)
set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)
xlabel('Time (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it u}-velocity (m/s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
h = legend(H, 'Flat fire', 'FDS', 'Location', 'Northeast');
set(h,'Interpreter', Font_Interpreter)
set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);
SVN_Filename = [repository, 'flat_fire_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end
display('Printing plot flat_fire_u.pdf...')
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/flat_fire_u')

close all
H(1) = plot(tvec, vexact, '-k');
hold on
H(2) = plot(STIME, QP(:,:,1,2), 'ok');
hold off
plot_style
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)
set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)
xlabel('Time (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it w}-velocity (m/s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
h = legend(H, 'Flat fire', 'FDS', 'Location', 'Southwest');
set(h,'Interpreter', Font_Interpreter)
set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);
SVN_Filename = [repository, 'flat_fire_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end
display('Printing plot flat_fire_w.pdf...')
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/flat_fire_w')

