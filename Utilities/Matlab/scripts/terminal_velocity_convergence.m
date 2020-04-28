% Trettel

close all
clear all

plot_style

repository = '../../Verification/Sprinklers_and_Sprays/';

skip_case = 0;

if ~exist([repository, 'terminal_velocity_dt_1_0_0001.prt5'])
    display(['Error: File ' [repository, 'terminal_velocity_dt_1_0_0001.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([repository, 'terminal_velocity_dt_0_1_0001.prt5'])
    display(['Error: File ' [repository, 'terminal_velocity_dt_0_1_0001.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([repository, 'terminal_velocity_dt_0_01_0001.prt5'])
    display(['Error: File ' [repository, 'terminal_velocity_dt_0_01_0001.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([repository, 'terminal_velocity_dt_0_001_0001.prt5'])
    display(['Error: File ' [repository, 'terminal_velocity_dt_0_001_0001.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([repository, 'terminal_velocity_dt_0_0001_0001.prt5'])
    display(['Error: File ' [repository, 'terminal_velocity_dt_0_0001_0001.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

% get precise ambient density
M = importdata([repository, 'terminal_velocity_dt_1_0_devc.csv'],',',2);
rhoa = M.data(end,2);

g = 9.8;
Cd = 1;
rhod = 1000;
D = 10e-3;

tend = 20;
eps = 1e-10;
ttest = 0.1;

K = 3 * rhoa * Cd / (4 * rhod * D);

errvec = [];
Linf = [];
vtexact = sqrt(g / K);
zexact = @(t) -log(cosh(sqrt(g * K) * t)) / K;
[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'terminal_velocity_dt_1_0_0001.prt5'],'real*8');
dtvec = [1.0 0.1 0.01 0.001 0.0001];
errvec(1) = abs(abs(QP(length(QP))) - vtexact);
Linf(1) = norm(ZP' - zexact(STIME), Inf);
%errtvec(1) = abs(ZP(find(abs(STIME - ttest) < eps, 1))' - zexact(ttest));
[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'terminal_velocity_dt_0_1_0001.prt5'],'real*8');
errvec(2) = abs(abs(QP(length(QP))) - vtexact);
Linf(2) = norm(ZP' - zexact(STIME), Inf);
%errtvec(2) = abs(ZP(find(abs(STIME - ttest) < eps, 1))' - zexact(ttest));
[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'terminal_velocity_dt_0_01_0001.prt5'],'real*8');
errvec(3) = abs(abs(QP(length(QP))) - vtexact);
Linf(3) = norm(ZP' - zexact(STIME), Inf);
%errtvec(3) = abs(ZP(find(abs(STIME - ttest) < eps, 1))' - zexact(ttest));
[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'terminal_velocity_dt_0_001_0001.prt5'],'real*8');
errvec(4) = abs(abs(QP(length(QP))) - vtexact);
Linf(4) = norm(ZP' - zexact(STIME), Inf);
%errtvec(4) = abs(ZP(find(abs(STIME - ttest) < eps, 1))' - zexact(ttest));
[STIME, XP, YP, ZP, QP] = read_prt5([repository, 'terminal_velocity_dt_0_0001_0001.prt5'],'real*8');
errvec(5) = abs(abs(QP(length(QP))) - vtexact);
Linf(5) = norm(ZP' - zexact(STIME), Inf);
%errtvec(6) = abs(ZP(find(abs(STIME - ttest) < eps, 1))' - zexact(ttest));

if errvec(5) > 1e-6
   display(['Matlab Warning: The velocity in the terminal_velocity* cases is out of tolerance.'])
end

if Linf(5) > 1e-6
   display(['Matlab Warning: The position in the terminal_velocity* cases is out of tolerance.'])
end

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1) = loglog(dtvec, errvec, '-*k');
hold on
H(2) = loglog(dtvec, dtvec, '--k');
H(3) = loglog(dtvec, dtvec.^2, '-k');
hold off
dto = dtvec((length(dtvec) - 1):length(dtvec));
erro = errvec((length(errvec) - 1):length(errvec));
[a, b] = polyfit(log(dto), log(erro), 1);
%fprintf('order of accuracy: %f\n', a(1))

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)

xlabel('Time Step (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Terminal Velocity Error (m/s)','FontSize',Label_Font_Size)
h = legend(H, 'FDS', 'O({\it \deltat})',...
    'O({\it \deltat}^2)', 'Location', 'Southeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

Git_Filename = [repository, 'terminal_velocity_dt_1_0_git.txt'];
addverstr(gca,Git_Filename,'loglog')

print(gcf, '-dpdf', '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/terminal_velocity_convergence');

close all
clear H
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1) = loglog(dtvec, Linf, '-*k');
hold on
H(2) = loglog(dtvec, 50 * dtvec, '--k');
H(3) = loglog(dtvec, 5 * dtvec.^2, '-k');
hold off
dto = dtvec((length(dtvec) - 1):length(dtvec));
erro = Linf((length(Linf) - 1):length(Linf));
[a, b] = polyfit(log(dto), log(erro), 1);
%fprintf('order of accuracy: %f\n', a(1))

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Label_Font_Size)

xlabel('Time Step (s)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Position Error (m)','FontSize',Label_Font_Size)
h = legend(H, 'FDS', 'O({\it \deltat})',...
    'O({\it \deltat}^2)', 'Location', 'Southeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

Git_Filename = [repository, 'terminal_velocity_dt_1_0_git.txt'];
addverstr(gca,Git_Filename,'loglog')

print(gcf, '-dpdf', '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/position_convergence');

