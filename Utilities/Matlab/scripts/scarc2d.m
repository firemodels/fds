% SKilian
% 4-8-10
% scard2d.m

close all
clear all

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

data_dir = '../../../Verification/Pressure_Solver/';
plot_dir = '../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

FFT_1MESH  = csvread([data_dir,'scarc2d_fft_1mesh_devc.csv'],2,0);
FFT_64MESH = csvread([data_dir,'scarc2d_fft_64mesh_devc.csv'],2,0);

SCARC_CG_64MESH    = csvread([data_dir,'scarc2d_cg_64mesh_devc.csv'],2,0);
SCARC_BICG_64MESH  = csvread([data_dir,'scarc2d_bicg_64mesh_devc.csv'],2,0);
SCARC_GMG_64MESH   = csvread([data_dir,'scarc2d_gmg_64mesh_devc.csv'],2,0);

time_fft_1mesh  = FFT_1MESH(:,1);  pres_fft_1mesh  = FFT_1MESH(:,5);
time_fft_64mesh = FFT_64MESH(:,1); pres_fft_64mesh = FFT_64MESH(:,5);

time_scarc_cg_64mesh    = SCARC_CG_64MESH(:,1)  ; pres_scarc_cg_64mesh   = SCARC_CG_64MESH(:,5);
time_scarc_bicg_64mesh  = SCARC_BICG_64MESH(:,1); pres_scarc_bicg_64mesh = SCARC_BICG_64MESH(:,5);
time_scarc_gmg_64mesh   = SCARC_GMG_64MESH(:,1) ; pres_scarc_gmg_64mesh  = SCARC_GMG_64MESH(:,5);

H(1)=plot(time_fft_1mesh ,pres_fft_1mesh ,'k-','LineWidth',Line_Width); hold on
H(2)=plot(time_fft_64mesh,pres_fft_64mesh,'k--','LineWidth',Line_Width); hold on

H(3)=plot(time_scarc_cg_64mesh(1:10:end) ,pres_scarc_cg_64mesh(1:10:end)   ,'ro','LineWidth',Line_Width); hold on
H(4)=plot(time_scarc_bicg_64mesh(1:10:end),pres_scarc_bicg_64mesh(1:10:end) ,'g^','LineWidth',Line_Width); hold on
H(5)=plot(time_scarc_gmg_64mesh(1:10:end) ,pres_scarc_gmg_64mesh(1:10:end)  ,'b+','LineWidth',Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter','LaTeX')
ylabel('Pressure (Pa)','FontSize',Title_Font_Size,'Interpreter','LaTeX')
axis([0 1 -2 2])
legend_handle=legend(H,'FFT 1 mesh','FFT 64 meshes','ScaRC-CG 64 meshes','ScaRC-BICG 64 meshes','ScaRC-GMG 64 meshes','Location','EastOutside');
legend boxoff
set(legend_handle,'FontSize',Title_Font_Size,'Interpreter','LaTeX')
set(legend_handle,'Position',[6.3 1.75 3.7 1.5])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% Create the PDF files
 
PDF_Paper_Width = 1.5*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'scarc2d'])




