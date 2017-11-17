% drjfloyd
% 9-3-2016
% fm_datacenter_veltest.m

% Generate velocity comparison plots for the FM datacenter flow mapping tests

expdir = '../../../exp/FM_FPRF_Datacenter/';
outdir = '../../../out/FM_FPRF_Datacenter/FDS_Output_Files/';

% High flow test
[exp_data] = csvread([expdir,'fm_datacenter_veltest_high.csv'],1);
[fds_data] = csvread([outdir,'FM_Datacenter_Veltest_High_devc.csv'],13);
n_fds_data = size(fds_data,1);

% compute average velocity
for i = 13:207
   fds_avg(i-12) = mean(fds_data(1:n_fds_data,i));
end

fds_u = fds_avg(1:65);
fds_v = fds_avg(66:130);
fds_w = fds_avg(131:195);
fds_u_rms = fds_data(n_fds_data,196:260);
fds_v_rms = fds_data(n_fds_data,261:325);
fds_w_rms = fds_data(n_fds_data,326:390);
fds_tot = (fds_u.^2+fds_v.^2+fds_w.^2).^0.5;
fds_tot_rms=(((fds_u.*fds_u_rms).^2+(fds_v.*fds_v_rms).^2+(fds_w.*fds_w_rms).^2)./fds_tot.^2).^0.5;

% set plot error lines

x_err = -10:1:10;
y_err = (0.1773^2 + (0.05*x_err).^2+(0.06*x_err).^2).^0.5;
y_err_p = x_err + 2*y_err;
y_err_m = x_err - 2*y_err;

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,1),fds_u,exp_data(:,4),'ro');
errorbar(exp_data(:,1),fds_u,fds_u_rms,'ro');
xlim([-2 2])
ylim([-2 2])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured U-Velocity (m/s)'];
ytitle = ['Predicted U-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

git_file = [outdir,'FM_Datacenter_Veltest_High_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_u'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,2),fds_v,exp_data(:,5),'ro');
errorbar(exp_data(:,2),fds_v,fds_v_rms,'ro');
xlim([-5 5])
ylim([-5 5])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured V-Velocity (m/s)'];
ytitle = ['Predicted V-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_v'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,3),fds_w,exp_data(:,6),'ro');
errorbar(exp_data(:,3),fds_w,fds_w_rms,'ro');
xlim([-1 3])
ylim([-1 3])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured W-Velocity (m/s)'];
ytitle = ['Predicted W-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_w'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,7),fds_tot,exp_data(:,8),'ro');
errorbar(exp_data(:,7),fds_tot,fds_tot_rms,'ro');
xlim([0 5])
ylim([0 5])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured Total Velocity (m/s)'];
ytitle = ['Predicted Total Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_vel'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

% Low flow test
[exp_data] = csvread([expdir,'fm_datacenter_veltest_low.csv'],1);
[fds_data] = csvread([outdir,'FM_Datacenter_Veltest_Low_devc.csv'],13);
n_fds_data = size(fds_data,1);

% compute average velocity
for i = 13:195
   fds_avg(i-12) = mean(fds_data(1:n_fds_data,i));
end

fds_u = fds_avg(1:61);
fds_v = fds_avg(62:122);
fds_w = fds_avg(123:183);
fds_u_rms = fds_data(n_fds_data,184:244);
fds_v_rms = fds_data(n_fds_data,245:305);
fds_w_rms = fds_data(n_fds_data,306:366);
fds_tot = (fds_u.^2+fds_v.^2+fds_w.^2).^0.5;
fds_tot_rms=(((fds_u.*fds_u_rms).^2+(fds_v.*fds_v_rms).^2+(fds_w.*fds_w_rms).^2)./fds_tot.^2).^0.5;

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,1),fds_u,exp_data(:,4),'ro');
errorbar(exp_data(:,1),fds_u,fds_u_rms,'ro');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured U-Velocity (m/s)'];
ytitle = ['Predicted U-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

git_file = [outdir,'FM_Datacenter_Veltest_Low_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_u'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,2),fds_v,exp_data(:,5),'ro');
errorbar(exp_data(:,2),fds_v,fds_v_rms,'ro');
xlim([-1.5 1.5])
ylim([-1.5 1.5])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured V-Velocity (m/s)'];
ytitle = ['Predicted V-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_v'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,3),fds_w,exp_data(:,6),'ro');
errorbar(exp_data(:,3),fds_w,fds_w_rms,'ro');
xlim([-0.4 0.8])
ylim([-0.4 0.8])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured W-Velocity (m/s)'];
ytitle = ['Predicted W-Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_w'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX

hX=plot(x_err,x_err,'k-',x_err,y_err_p,'k--',x_err,y_err_m,'k--');
hold on
herrorbar(exp_data(:,7),fds_tot,exp_data(:,8),'ro');
errorbar(exp_data(:,7),fds_tot,fds_tot_rms,'ro');
xlim([0 1.4])
ylim([0 1.4])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured Total Velocity (m/s)'];
ytitle = ['Predicted Total Velocity (m/s)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_vel'];
warning('off','MATLAB:print:FigureTooLargeForPage')
print(gcf,'-dpdf',plotname);
hold off
clear hX


