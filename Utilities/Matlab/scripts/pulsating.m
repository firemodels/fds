% Antonellis
% pulsating.m
% 06-24-10

close all
clear all

plot_style

data_dir = '../../Verification/Scalar_Analytical_Solution/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

devc_col = 3; % (x,y)=(pi,pi)-->devc_col=2, (x,y)=(1.5*pi,1.5*pi)-->devc_col=3

skip_case = 0;

if ~exist([data_dir,'pulsating_FL2_16_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL2_32_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL2_64_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL2_128_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pulsating_FL4_16_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL4_32_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL4_64_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL4_128_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pulsating_FL0_16_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL0_32_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL0_64_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_64_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'pulsating_FL0_128_devc.csv'])
    display(['Error: File ' [data_dir,'pulsating_FL2_128_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

% Superbee limiter, FL=2
M_FL2_16 = csvread([data_dir,'pulsating_FL2_16_devc.csv'],2,0);
M_FL2_32 = csvread([data_dir,'pulsating_FL2_32_devc.csv'],2,0);
M_FL2_64 = csvread([data_dir,'pulsating_FL2_64_devc.csv'],2,0);
M_FL2_128 = csvread([data_dir,'pulsating_FL2_128_devc.csv'],2,0);
t_FL2_16 = M_FL2_16(:,1); rho_fds_FL2_16 = M_FL2_16(:,devc_col);
t_FL2_32 = M_FL2_32(:,1); rho_fds_FL2_32 = M_FL2_32(:,devc_col);
t_FL2_64 = M_FL2_64(:,1); rho_fds_FL2_64 = M_FL2_64(:,devc_col);
t_FL2_128 = M_FL2_128(:,1); rho_fds_FL2_128 = M_FL2_128(:,devc_col);

% CHARM limiter, FL=4
M_FL4_16 = csvread([data_dir,'pulsating_FL4_16_devc.csv'],2,0);
M_FL4_32 = csvread([data_dir,'pulsating_FL4_32_devc.csv'],2,0);
M_FL4_64 = csvread([data_dir,'pulsating_FL4_64_devc.csv'],2,0);
M_FL4_128 = csvread([data_dir,'pulsating_FL4_128_devc.csv'],2,0);
t_FL4_16 = M_FL4_16(:,1); rho_fds_FL4_16 = M_FL4_16(:,devc_col);
t_FL4_32 = M_FL4_32(:,1); rho_fds_FL4_32 = M_FL4_32(:,devc_col);
t_FL4_64 = M_FL4_64(:,1); rho_fds_FL4_64 = M_FL4_64(:,devc_col);
t_FL4_128 = M_FL4_128(:,1); rho_fds_FL4_128 = M_FL4_128(:,devc_col);

% central differencing, FL=0
M_FL0_16 = csvread([data_dir,'pulsating_FL0_16_devc.csv'],2,0);
M_FL0_32 = csvread([data_dir,'pulsating_FL0_32_devc.csv'],2,0);
M_FL0_64 = csvread([data_dir,'pulsating_FL0_64_devc.csv'],2,0);
M_FL0_128 = csvread([data_dir,'pulsating_FL0_128_devc.csv'],2,0);
t_FL0_16 = M_FL0_16(:,1); rho_fds_FL0_16 = M_FL0_16(:,devc_col);
t_FL0_32 = M_FL0_32(:,1); rho_fds_FL0_32 = M_FL0_32(:,devc_col);
t_FL0_64 = M_FL0_64(:,1); rho_fds_FL0_64 = M_FL0_64(:,devc_col);
t_FL0_128 = M_FL0_128(:,1); rho_fds_FL0_128 = M_FL0_128(:,devc_col);

% analytical solution

B = 1;
w = 1;

L = 2*pi;
if devc_col==2; x=pi; y=pi; end
if devc_col==3; x=1.5*pi; y=1.5*pi; end

rho_FL2_16 = section2_soln(rho_fds_FL2_16(1),x-L/32,y-L/32,B,w,t_FL2_16);      error_FL2_16 = norm(rho_fds_FL2_16-rho_FL2_16)/length(t_FL2_16);
rho_FL2_32 = section2_soln(rho_fds_FL2_32(1),x-L/64,y-L/64,B,w,t_FL2_32);      error_FL2_32 = norm(rho_fds_FL2_32-rho_FL2_32)/length(t_FL2_32);
rho_FL2_64 = section2_soln(rho_fds_FL2_64(1),x-L/128,y-L/128,B,w,t_FL2_64);    error_FL2_64 = norm(rho_fds_FL2_64-rho_FL2_64)/length(t_FL2_64);
rho_FL2_128 = section2_soln(rho_fds_FL2_128(1),x-L/256,y-L/256,B,w,t_FL2_128); error_FL2_128 = norm(rho_fds_FL2_128-rho_FL2_128)/length(t_FL2_128);

rho_FL4_16 = section2_soln(rho_fds_FL4_16(1),x-L/32,y-L/32,B,w,t_FL4_16);      error_FL4_16 = norm(rho_fds_FL4_16-rho_FL4_16)/length(t_FL4_16);
rho_FL4_32 = section2_soln(rho_fds_FL4_32(1),x-L/64,y-L/64,B,w,t_FL4_32);      error_FL4_32 = norm(rho_fds_FL4_32-rho_FL4_32)/length(t_FL4_32);
rho_FL4_64 = section2_soln(rho_fds_FL4_64(1),x-L/128,y-L/128,B,w,t_FL4_64);    error_FL4_64 = norm(rho_fds_FL4_64-rho_FL4_64)/length(t_FL4_64);
rho_FL4_128 = section2_soln(rho_fds_FL4_128(1),x-L/256,y-L/256,B,w,t_FL4_128); error_FL4_128 = norm(rho_fds_FL4_128-rho_FL4_128)/length(t_FL4_128);

rho_FL0_16 = section2_soln(rho_fds_FL0_16(1),x-L/32,y-L/32,B,w,t_FL0_16);      error_FL0_16 = norm(rho_fds_FL0_16-rho_FL0_16)/length(t_FL0_16);
rho_FL0_32 = section2_soln(rho_fds_FL0_32(1),x-L/64,y-L/64,B,w,t_FL0_32);      error_FL0_32 = norm(rho_fds_FL0_32-rho_FL0_32)/length(t_FL0_32);
rho_FL0_64 = section2_soln(rho_fds_FL0_64(1),x-L/128,y-L/128,B,w,t_FL0_64);    error_FL0_64 = norm(rho_fds_FL0_64-rho_FL0_64)/length(t_FL0_64);
rho_FL0_128 = section2_soln(rho_fds_FL0_128(1),x-L/256,y-L/256,B,w,t_FL0_128); error_FL0_128 = norm(rho_fds_FL0_128-rho_FL0_128)/length(t_FL0_128);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(5)=plot(t_FL2_128,rho_FL2_128,'k-','LineWidth',Line_Width); hold on
H(1)=plot(t_FL2_16,rho_fds_FL2_16,'c--','LineWidth',Line_Width); hold on
H(2)=plot(t_FL2_32,rho_fds_FL2_32,'g--','LineWidth',Line_Width); hold on
H(3)=plot(t_FL2_64,rho_fds_FL2_64,'b--','LineWidth',Line_Width); hold on
H(4)=plot(t_FL2_128,rho_fds_FL2_128,'r--','LineWidth',Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('Density (kg/m^3)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([0 12.5 0 2])
legend_handle=legend(H,'FDS N=16','FDS N=32','FDS N=64','FDS N=128','Analytical Solution','Location','NorthEast');
set(legend_handle,'Interpreter',Font_Interpreter);
set(legend_handle,'Fontsize',Key_Font_Size);
set(legend_handle,'Box','on');
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% add Git revision if file is available

Git_Filename = [data_dir,'pulsating_FL2_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pulsating_time_series'])

% convergence plot

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h = 2*pi./[16 32 64 128];
e_FL2 = [error_FL2_16 error_FL2_32 error_FL2_64 error_FL2_128];
e_FL4 = [error_FL4_16 error_FL4_32 error_FL4_64 error_FL4_128];
e_FL0 = [error_FL0_16 error_FL0_32 error_FL0_64 error_FL0_128];

H(1)=loglog(h,e_FL0,'b*-','LineWidth',Line_Width); hold on
H(2)=loglog(h,e_FL2,'ro-','LineWidth',Line_Width); hold on
H(3)=loglog(h,e_FL4,'g^-','LineWidth',Line_Width); hold on
H(4)=loglog(h,.1*h,'k--','LineWidth',Line_Width);
H(5)=loglog(h,.1*h.^2,'k-','LineWidth',Line_Width);

xlabel('Grid Spacing (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('L2 Error (kg/m^3)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([1e-2 1e0 1e-6 1e-1])
legend_handle=legend(H(1:5),'FDS Central','FDS Superbee','FDS CHARM','{\it O}({\it\deltax})','{\it O}({\it\deltax^2})','Location','NorthWest');
set(legend_handle,'Interpreter',Font_Interpreter);
set(legend_handle,'Fontsize',Key_Font_Size);
set(legend_handle,'Box','on');
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% add Git revision if file is available

Git_Filename = [data_dir,'pulsating_FL2_16_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pulsating_convergence'])



