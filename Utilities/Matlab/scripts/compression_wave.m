% McDermott
% 4-8-10
% compression_wave.m

close all
clear all

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

data_dir = '../../Verification/Scalar_Analytical_Solution/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

skip_case = 0;

if ~exist([data_dir,'compression_wave_FL0_16_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL0_32_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL0_64_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL0_128_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'compression_wave_FL2_16_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL2_32_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL2_64_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL2_128_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'compression_wave_FL4_16_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL4_32_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL4_64_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL4_128_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'compression_wave_FL5_16_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL5_32_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL5_64_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'compression_wave_FL5_128_devc.csv'])
    display(['Error: File ' [data_dir,'compression_wave_FL0_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

% central differencing, FL=0

M_FL0_16 = csvread([data_dir,'compression_wave_FL0_16_devc.csv'],2,0);
M_FL0_32 = csvread([data_dir,'compression_wave_FL0_32_devc.csv'],2,0);
M_FL0_64 = csvread([data_dir,'compression_wave_FL0_64_devc.csv'],2,0);
M_FL0_128 = csvread([data_dir,'compression_wave_FL0_128_devc.csv'],2,0);
t_FL0_16 = M_FL0_16(:,1); rho_fds_FL0_16 = M_FL0_16(:,3);
t_FL0_32 = M_FL0_32(:,1); rho_fds_FL0_32 = M_FL0_32(:,3);
t_FL0_64 = M_FL0_64(:,1); rho_fds_FL0_64 = M_FL0_64(:,3);
t_FL0_128 = M_FL0_128(:,1); rho_fds_FL0_128 = M_FL0_128(:,3);

% Superbee limiter, FL=2
M_FL2_16 = csvread([data_dir,'compression_wave_FL2_16_devc.csv'],2,0);
M_FL2_32 = csvread([data_dir,'compression_wave_FL2_32_devc.csv'],2,0);
M_FL2_64 = csvread([data_dir,'compression_wave_FL2_64_devc.csv'],2,0);
M_FL2_128 = csvread([data_dir,'compression_wave_FL2_128_devc.csv'],2,0);
t_FL2_16 = M_FL2_16(:,1); rho_fds_FL2_16 = M_FL2_16(:,3);
t_FL2_32 = M_FL2_32(:,1); rho_fds_FL2_32 = M_FL2_32(:,3);
t_FL2_64 = M_FL2_64(:,1); rho_fds_FL2_64 = M_FL2_64(:,3);
t_FL2_128 = M_FL2_128(:,1); rho_fds_FL2_128 = M_FL2_128(:,3);

% CHARM limiter, FL=4
M_FL4_16 = csvread([data_dir,'compression_wave_FL4_16_devc.csv'],2,0);
M_FL4_32 = csvread([data_dir,'compression_wave_FL4_32_devc.csv'],2,0);
M_FL4_64 = csvread([data_dir,'compression_wave_FL4_64_devc.csv'],2,0);
M_FL4_128 = csvread([data_dir,'compression_wave_FL4_128_devc.csv'],2,0);
t_FL4_16 = M_FL4_16(:,1); rho_fds_FL4_16 = M_FL4_16(:,3);
t_FL4_32 = M_FL4_32(:,1); rho_fds_FL4_32 = M_FL4_32(:,3);
t_FL4_64 = M_FL4_64(:,1); rho_fds_FL4_64 = M_FL4_64(:,3);
t_FL4_128 = M_FL4_128(:,1); rho_fds_FL4_128 = M_FL4_128(:,3);

% MP5 limiter, FL=5
M_FL5_16 = csvread([data_dir,'compression_wave_FL5_16_devc.csv'],2,0);
M_FL5_32 = csvread([data_dir,'compression_wave_FL5_32_devc.csv'],2,0);
M_FL5_64 = csvread([data_dir,'compression_wave_FL5_64_devc.csv'],2,0);
M_FL5_128 = csvread([data_dir,'compression_wave_FL5_128_devc.csv'],2,0);
t_FL5_16 = M_FL5_16(:,1); rho_fds_FL5_16 = M_FL5_16(:,3);
t_FL5_32 = M_FL5_32(:,1); rho_fds_FL5_32 = M_FL5_32(:,3);
t_FL5_64 = M_FL5_64(:,1); rho_fds_FL5_64 = M_FL5_64(:,3);
t_FL5_128 = M_FL5_128(:,1); rho_fds_FL5_128 = M_FL5_128(:,3);

% analytical solution

a = 2;
c = 3;

L = 2*pi;
x = 1.5*pi;
y = 1.5*pi;

rho_FL0_16 = compression_wave_soln(rho_fds_FL0_16(1),x-L/32,y-L/32,a,c,t_FL0_16);      error_FL0_16 = norm(rho_fds_FL0_16-rho_FL0_16)/length(t_FL0_16);
rho_FL0_32 = compression_wave_soln(rho_fds_FL0_32(1),x-L/64,y-L/64,a,c,t_FL0_32);      error_FL0_32 = norm(rho_fds_FL0_32-rho_FL0_32)/length(t_FL0_32);
rho_FL0_64 = compression_wave_soln(rho_fds_FL0_64(1),x-L/128,y-L/128,a,c,t_FL0_64);    error_FL0_64 = norm(rho_fds_FL0_64-rho_FL0_64)/length(t_FL0_64);
rho_FL0_128 = compression_wave_soln(rho_fds_FL0_128(1),x-L/256,y-L/256,a,c,t_FL0_128); error_FL0_128 = norm(rho_fds_FL0_128-rho_FL0_128)/length(t_FL0_128);

rho_FL2_16 = compression_wave_soln(rho_fds_FL2_16(1),x-L/32,y-L/32,a,c,t_FL2_16);      error_FL2_16 = norm(rho_fds_FL2_16-rho_FL2_16)/length(t_FL2_16);
rho_FL2_32 = compression_wave_soln(rho_fds_FL2_32(1),x-L/64,y-L/64,a,c,t_FL2_32);      error_FL2_32 = norm(rho_fds_FL2_32-rho_FL2_32)/length(t_FL2_32);
rho_FL2_64 = compression_wave_soln(rho_fds_FL2_64(1),x-L/128,y-L/128,a,c,t_FL2_64);    error_FL2_64 = norm(rho_fds_FL2_64-rho_FL2_64)/length(t_FL2_64);
rho_FL2_128 = compression_wave_soln(rho_fds_FL2_128(1),x-L/256,y-L/256,a,c,t_FL2_128); error_FL2_128 = norm(rho_fds_FL2_128-rho_FL2_128)/length(t_FL2_128);

rho_FL4_16 = compression_wave_soln(rho_fds_FL4_16(1),x-L/32,y-L/32,a,c,t_FL4_16);      error_FL4_16 = norm(rho_fds_FL4_16-rho_FL4_16)/length(t_FL4_16);
rho_FL4_32 = compression_wave_soln(rho_fds_FL4_32(1),x-L/64,y-L/64,a,c,t_FL4_32);      error_FL4_32 = norm(rho_fds_FL4_32-rho_FL4_32)/length(t_FL4_32);
rho_FL4_64 = compression_wave_soln(rho_fds_FL4_64(1),x-L/128,y-L/128,a,c,t_FL4_64);    error_FL4_64 = norm(rho_fds_FL4_64-rho_FL4_64)/length(t_FL4_64);
rho_FL4_128 = compression_wave_soln(rho_fds_FL4_128(1),x-L/256,y-L/256,a,c,t_FL4_128); error_FL4_128 = norm(rho_fds_FL4_128-rho_FL4_128)/length(t_FL4_128);

rho_FL5_16 = compression_wave_soln(rho_fds_FL5_16(1),x-L/32,y-L/32,a,c,t_FL5_16);      error_FL5_16 = norm(rho_fds_FL5_16-rho_FL5_16)/length(t_FL5_16);
rho_FL5_32 = compression_wave_soln(rho_fds_FL5_32(1),x-L/64,y-L/64,a,c,t_FL5_32);      error_FL5_32 = norm(rho_fds_FL5_32-rho_FL5_32)/length(t_FL5_32);
rho_FL5_64 = compression_wave_soln(rho_fds_FL5_64(1),x-L/128,y-L/128,a,c,t_FL5_64);    error_FL5_64 = norm(rho_fds_FL5_64-rho_FL5_64)/length(t_FL5_64);
rho_FL5_128 = compression_wave_soln(rho_fds_FL5_128(1),x-L/256,y-L/256,a,c,t_FL5_128); error_FL5_128 = norm(rho_fds_FL5_128-rho_FL5_128)/length(t_FL5_128);


H(5)=plot(t_FL4_128,rho_FL4_128,'k-','LineWidth',Line_Width); hold on
H(1)=plot(t_FL4_16,rho_fds_FL4_16,'c--','LineWidth',Line_Width); hold on
H(2)=plot(t_FL4_32,rho_fds_FL4_32,'g--','LineWidth',Line_Width); hold on
H(3)=plot(t_FL4_64,rho_fds_FL4_64,'b--','LineWidth',Line_Width); hold on
H(4)=plot(t_FL4_128,rho_fds_FL4_128,'r--','LineWidth',Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('Density (kg/m^3)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([0 12.5 0 8])
legend_handle=legend(H,'FDS N=16','FDS N=32','FDS N=64','FDS N=128','Analytical Solution','Location','NorthEast');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% add SVN if file is available

SVN_Filename = [data_dir,'compression_wave_FL0_16_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'compression_wave_time_series'])

% convergence plot

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

h = 2*pi./[16 32 64 128];
e_FL0 = [error_FL0_16 error_FL0_32 error_FL0_64 error_FL0_128];
e_FL2 = [error_FL2_16 error_FL2_32 error_FL2_64 error_FL2_128];
e_FL4 = [error_FL4_16 error_FL4_32 error_FL4_64 error_FL4_128];
e_FL5 = [error_FL5_16 error_FL5_32 error_FL5_64 error_FL5_128];
H(1)=loglog(h,e_FL0,'b*-','LineWidth',Line_Width); hold on
H(2)=loglog(h,e_FL2,'ro-','LineWidth',Line_Width); hold on
H(3)=loglog(h,e_FL4,'g^-','LineWidth',Line_Width); hold on
H(4)=loglog(h,e_FL5,'c>-','LineWidth',Line_Width); hold on
H(5)=loglog(h,.1*h,'k--','LineWidth',Line_Width);
H(6)=loglog(h,.1*h.^2,'k-','LineWidth',Line_Width);

xlabel('Grid Spacing (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('L2 Error (kg/m^3)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([1e-2 1e0 1e-4 1e-1])
legend_handle=legend(H(1:6),'Central','Superbee','CHARM','MP5','{\it O}({\it\deltax})','{\it O}({\it\deltax^2})','Location','NorthWest');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% add SVN if file is available

SVN_Filename = [data_dir,'compression_wave_FL0_16_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'compression_wave_convergence'])



