% McDermott
% 1-9-2017
% umd_line_burner_2.m
%
% Plot combustion efficiency for UMD Line Burner cases.

close all
clear all

plot_style

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';


EXP = importdata([expdir,'XO2_vs_ETA_METHANE.csv'],',',1);
%HRR = importdata([outdir,'methane_dx_XO2_ramp_dx_p625cm_hrr.csv'],',',2);

XO2 = EXP.data(:,find(strcmp(EXP.colheaders,'XO2')));
eta = EXP.data(:,find(strcmp(EXP.colheaders,'eta')));

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(XO2,eta,'k.','MarkerSize',10);
axis([0.09  0.21 0 1.1 ])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\eta','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exp','Location','NorthWest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_eta']);