% SKilian
% 4-8-10
% scard2d.m

close all
clear all

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

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

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Time (s)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Pressure (Pa)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
axis([0 1 -2 2])

legend_handle=legend(H,'FFT 1 mesh','FFT 8 meshes','ScaRC-CG 8 meshes','ScaRC-BICG 8 meshes','ScaRC-GMG 8 meshes','Location','EastOutside');
set(legend_handle,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)
set(legend_handle,'Units',Paper_Units)
pos = get(legend_handle,'position');
set(legend_handle,'position',[Paper_Width pos(2:4)])

Git_Filename = [data_dir,'scarc2d_cg_8mesh_git.txt'];
addverstr(gca,Git_Filename,'linear')

% Create the PDF files

PDF_Paper_Width = 1.55*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'scarc2d'])




