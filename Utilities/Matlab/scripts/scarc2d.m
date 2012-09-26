% SKilian
% 4-8-10
% scard2d.m

close all
clear all

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

data_dir = '../../Verification/Pressure_Solver/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

skip_case = 0;

if ~exist([data_dir,'scarc2d_fft_1mesh_devc.csv'])
    display(['Error: File ' [data_dir,'scarc2d_fft_1mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'scarc2d_fft_8mesh_devc.csv'])
    display(['Error: File ' [data_dir,'scarc2d_fft_8mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'scarc2d_cg_8mesh_devc.csv'])
    display(['Error: File ' [data_dir,'scarc2d_cg_8mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'scarc2d_bicg_8mesh_devc.csv'])
    display(['Error: File ' [data_dir,'scarc2d_bicg_8mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([data_dir,'scarc2d_gmg_8mesh_devc.csv'])
    display(['Error: File ' [data_dir,'scarc2d_gmg_8mesh_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

FFT_1MESH  = csvread([data_dir,'scarc2d_fft_1mesh_devc.csv'],2,0);
FFT_8MESH = csvread([data_dir,'scarc2d_fft_8mesh_devc.csv'],2,0);

SCARC_CG_8MESH    = csvread([data_dir,'scarc2d_cg_8mesh_devc.csv'],2,0);
SCARC_BICG_8MESH  = csvread([data_dir,'scarc2d_bicg_8mesh_devc.csv'],2,0);
SCARC_GMG_8MESH   = csvread([data_dir,'scarc2d_gmg_8mesh_devc.csv'],2,0);

time_fft_1mesh  = FFT_1MESH(:,1);  pres_fft_1mesh  = FFT_1MESH(:,5);
time_fft_8mesh = FFT_8MESH(:,1); pres_fft_8mesh = FFT_8MESH(:,5);

time_scarc_cg_8mesh    = SCARC_CG_8MESH(:,1)  ; pres_scarc_cg_8mesh   = SCARC_CG_8MESH(:,5);
time_scarc_bicg_8mesh  = SCARC_BICG_8MESH(:,1); pres_scarc_bicg_8mesh = SCARC_BICG_8MESH(:,5);
time_scarc_gmg_8mesh   = SCARC_GMG_8MESH(:,1) ; pres_scarc_gmg_8mesh  = SCARC_GMG_8MESH(:,5);

H(1)=plot(time_fft_1mesh ,pres_fft_1mesh ,'k-','LineWidth',Line_Width); hold on
H(2)=plot(time_fft_8mesh,pres_fft_8mesh,'k--','LineWidth',Line_Width); hold on

H(3)=plot(time_scarc_cg_8mesh(1:10:end) ,pres_scarc_cg_8mesh(1:10:end)   ,'ro','LineWidth',Line_Width); hold on
H(4)=plot(time_scarc_bicg_8mesh(1:10:end),pres_scarc_bicg_8mesh(1:10:end) ,'g^','LineWidth',Line_Width); hold on
H(5)=plot(time_scarc_gmg_8mesh(1:10:end) ,pres_scarc_gmg_8mesh(1:10:end)  ,'b+','LineWidth',Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Pressure (Pa)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
axis([0 1 -2 2])
legend_handle=legend(H,'FFT 1 mesh','FFT 8 meshes','ScaRC-CG 8 meshes','ScaRC-BICG 8 meshes','ScaRC-GMG 8 meshes','Location','EastOutside');
legend boxoff
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
set(legend_handle,'Position',[6.3 1.75 3.7 1.5])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

SVN_Filename = [data_dir,'scarc2d_cg_8mesh_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% Create the PDF files
 
PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'scarc2d'])




