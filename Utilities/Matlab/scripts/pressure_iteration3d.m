% SKilian
% 3-20-18
% pressure_iteration3d.m

close all
clear all

data_dir = '../../../Verification/Pressure_Solver/';
plot_dir = '../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

SVN_Filename = [data_dir,'pressure_iteration3d_8mesh_git.txt'];

skip_case = 0;

if ~exist([data_dir,'pressure_iteration3d_1mesh_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_1mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_default_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_default_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_tight_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_tight_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_glmat_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_glmat_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_uscarc_cg_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_uscarc_cg_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_scarc_cg_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_scarc_cg_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_scarc_cg_tight_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_scarc_cg_tight_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_scarc_gmg_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_scarc_gmg_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_scarc_cggmg_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_scarc_cggmg_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([data_dir,'pressure_iteration3d_scarc_twolevel_devc.csv'])
    display(['Error: File ' [data_dir,'pressure_iteration3d_scarc_twolevel.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

FFT_1mesh      = csvread([data_dir,'pressure_iteration3d_1mesh_devc.csv'],2,0);
FFT_default    = csvread([data_dir,'pressure_iteration3d_default_devc.csv'],2,0);
FFT_tight      = csvread([data_dir,'pressure_iteration3d_tight_devc.csv'],2,0);
GLMAT          = csvread([data_dir,'pressure_iteration3d_glmat_devc.csv'],2,0);
USCARC_cg      = csvread([data_dir,'pressure_iteration3d_uscarc_cg_devc.csv'],2,0);
SCARC_cg       = csvread([data_dir,'pressure_iteration3d_scarc_cg_devc.csv'],2,0);
SCARC_cg_tight = csvread([data_dir,'pressure_iteration3d_scarc_cg_tight_devc.csv'],2,0);
SCARC_gmg      = csvread([data_dir,'pressure_iteration3d_scarc_gmg_devc.csv'],2,0);
SCARC_cggmg    = csvread([data_dir,'pressure_iteration3d_scarc_cggmg_devc.csv'],2,0);
SCARC_twolevel = csvread([data_dir,'pressure_iteration3d_scarc_twolevel_devc.csv'],2,0);

time_fft_1mesh      = FFT_1mesh(:,1);       
time_fft_default    = FFT_default(:,1);    
time_fft_tight      = FFT_tight(:,1);      
time_glmat          = GLMAT(:,1);          
time_uscarc_cg      = USCARC_cg(:,1);      
time_scarc_cg       = SCARC_cg(:,1);       
time_scarc_cg_tight = SCARC_cg_tight(:,1); 
time_scarc_gmg      = SCARC_gmg(:,1);      
time_scarc_cggmg    = SCARC_cggmg(:,1);    
time_scarc_twolevel = SCARC_twolevel(:,1); 


% pressure plot 

figure(1)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

pres_fft_1mesh      = FFT_1mesh(:,5);       
pres_fft_default    = FFT_default(:,5);    
pres_fft_tight      = FFT_tight(:,5);      
pres_glmat          = GLMAT(:,5);          
pres_scarc_cg       = SCARC_cg(:,5);       
pres_uscarc_cg      = USCARC_cg(:,5);      

H(1) = plot(time_fft_1mesh,   pres_fft_1mesh,   'r',   'LineWidth', Line_Width); hold on
H(2) = plot(time_fft_default, pres_fft_default, 'k',   'LineWidth', Line_Width); hold on
H(3) = plot(time_fft_tight,   pres_fft_tight,   'c',   'LineWidth', Line_Width); hold on
H(4) = plot(time_glmat,       pres_glmat,       'm',   'LineWidth', Line_Width); hold on
H(5) = plot(time_scarc_cg,    pres_scarc_cg,    'b-.', 'LineWidth', Line_Width); hold on
H(6) = plot(time_uscarc_cg,   pres_uscarc_cg,   'g--', 'LineWidth', Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Pressure (Pa)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
axis([0 0.5 -0.4 0.35])
legend_handle=legend(H,'FFT 1 mesh','FFT default','FFT tight','GLMAT','ScaRC','UScaRC','Location','EastOutside');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

addverstr(gca,SVN_Filename,'linear')

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pressure_iteration3d_pres'])



% velocity error plot

figure(2)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

verr_fft_1mesh      = FFT_1mesh(:,6);       
verr_fft_default    = FFT_default(:,6);    
verr_fft_tight      = FFT_tight(:,6);      
verr_glmat          = GLMAT(:,6);          
verr_uscarc_cg      = USCARC_cg(:,6);      
verr_scarc_cg       = SCARC_cg(:,6);       
verr_scarc_cg_tight = SCARC_cg_tight(:,6); 

H(1) = plot(time_fft_1mesh,   verr_fft_1mesh,   'r',   'LineWidth', Line_Width); hold on
H(2) = plot(time_fft_default, verr_fft_default, 'k',   'LineWidth', Line_Width); hold on
H(3) = plot(time_fft_tight,   verr_fft_tight,   'c',   'LineWidth', Line_Width); hold on
H(4) = plot(time_glmat,       verr_glmat,       'm',   'LineWidth', Line_Width); hold on
H(5) = plot(time_scarc_cg,    verr_scarc_cg,    'b-.', 'LineWidth', Line_Width); hold on
H(6) = plot(time_uscarc_cg,   verr_uscarc_cg,   'g--', 'LineWidth', Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Velocity error (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
%text(0.1, 0.1, 'Velocity error ','FontSize',Label_Font_Size,'FontName',Font_Name)
axis([0 0.5 1e-17 1])
legend_handle=legend(H,'FFT 1 mesh','FFT default','FFT tight','GLMAT', 'ScaRC','UScaRC','Location','EastOutside');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);
set(gca, 'YScale', 'log')
set(gca,'FontSize',Title_Font_Size)

addverstr(gca,SVN_Filename,'linear')

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pressure_iteration3d_verr'])


% pressure iterations plot

figure(3)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

pite_fft_1mesh      = FFT_1mesh(:,7);       
pite_fft_default    = FFT_default(:,7);    
pite_fft_tight      = FFT_tight(:,7);      
pite_glmat          = GLMAT(:,7);          
pite_uscarc_cg      = USCARC_cg(:,7);      
pite_scarc_cg       = SCARC_cg(:,7);       
pite_scarc_cg_tight = SCARC_cg_tight(:,7); 

H(1) = plot(time_fft_1mesh,   pite_fft_1mesh,   'r',   'LineWidth', Line_Width); hold on
H(2) = plot(time_fft_default, pite_fft_default, 'k',   'LineWidth', Line_Width); hold on
H(3) = plot(time_fft_tight,   pite_fft_tight,   'c',   'LineWidth', Line_Width); hold on
H(4) = plot(time_glmat,       pite_glmat,       'm',   'LineWidth', Line_Width); hold on
H(5) = plot(time_scarc_cg,    pite_scarc_cg,    'b-.', 'LineWidth', Line_Width); hold on
H(6) = plot(time_uscarc_cg,   pite_uscarc_cg,   'g--', 'LineWidth', Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Pressure iterations','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
axis([0 0.5 0 11])
legend_handle=legend(H,'FFT 1 mesh','FFT default','FFT tight','GLMAT', 'ScaRC','UScaRC','Location','EastOutside');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

addverstr(gca,SVN_Filename,'linear')

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pressure_iteration3d_pite'])


% ScaRC iterations plot

figure(4)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

site_uscarc_cg      = USCARC_cg(:,8);      
site_scarc_cg       = SCARC_cg(:,8);       
site_scarc_cg_tight = SCARC_cg_tight(:,8); 
site_scarc_gmg      = SCARC_gmg(:,8);      
site_scarc_cggmg    = SCARC_cggmg(:,8);    
site_scarc_twolevel = SCARC_twolevel(:,8); 

H(1) = plot(time_uscarc_cg,      site_uscarc_cg,      'k', 'LineWidth', Line_Width); hold on
H(2) = plot(time_scarc_cg,       site_scarc_cg,       'b', 'LineWidth', Line_Width); hold on
H(3) = plot(time_scarc_cg_tight, site_scarc_cg_tight, 'm', 'LineWidth', Line_Width); hold on
H(4) = plot(time_scarc_gmg,      site_scarc_gmg,      'c', 'LineWidth', Line_Width); hold on
H(5) = plot(time_scarc_cggmg,    site_scarc_cggmg,    'r', 'LineWidth', Line_Width); hold on
H(6) = plot(time_scarc_twolevel, site_scarc_twolevel, 'g', 'LineWidth', Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('ScaRC iterations','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
axis([0 0.5 0 55])
legend_handle=legend(H,'UScaRC','ScaRC-CG','ScaRC-CG tight','ScaRC-GMG','ScaRC-CGGMG','ScaRC-Twolevel', 'Location','EastOutside');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

addverstr(gca,SVN_Filename,'linear')

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pressure_iteration3d_site'])


% ScaRC convergence rate plot

figure(5)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

scon_uscarc_cg      = USCARC_cg(:,9);      
scon_scarc_cg       = SCARC_cg(:,9);       
scon_scarc_cg_tight = SCARC_cg_tight(:,9); 
scon_scarc_gmg      = SCARC_gmg(:,9);      
scon_scarc_cggmg    = SCARC_cggmg(:,9);    
scon_scarc_twolevel = SCARC_twolevel(:,9); 

H(1) = plot(time_uscarc_cg,      scon_uscarc_cg,      'k', 'LineWidth', Line_Width); hold on
H(2) = plot(time_scarc_cg,       scon_scarc_cg,       'b', 'LineWidth', Line_Width); hold on
H(3) = plot(time_scarc_cg_tight, scon_scarc_cg_tight, 'm', 'LineWidth', Line_Width); hold on
H(4) = plot(time_scarc_gmg,      scon_scarc_gmg,      'c', 'LineWidth', Line_Width); hold on
H(5) = plot(time_scarc_cggmg,    scon_scarc_cggmg,    'r', 'LineWidth', Line_Width); hold on
H(6) = plot(time_scarc_twolevel, scon_scarc_twolevel, 'g', 'LineWidth', Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('ScaRC convergence rate','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
axis([0 0.5 0 0.7])
legend_handle=legend(H,'UScaRC', 'ScaRC-CG','ScaRC-CG tight','ScaRC-GMG','ScaRC-CGGMG','ScaRC-Twolevel','Location','EastOutside');
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

addverstr(gca,SVN_Filename,'linear')

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'pressure_iteration3d_scon'])



