% Topi Sikanen
% Oct 25, 2013
% VTT_Sprays.m

disp('VTT_Sprays ...')

plot_style

outdir = '../../../out/VTT_Sprays/';
expdir = '../../../exp/VTT_Sprays/';

files = dir(strcat([outdir,'*_line.csv']));

skip_case = 0;
if ~exist([outdir,'LN02_1_line.csv'])
    display('Error: File LN02_1_line.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'LN02_2_line.csv'])
    display('Error: File LN02_2_line.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'LN02_4_line.csv'])
    display('Error: File LN02_4_line.csv does not exist. Skipping case.')
    skip_case = 1;
end
if skip_case
    return
end

edata=importdata([expdir,'LN02.csv']);
L1=importdata([outdir,'LN02_1_line.csv']);
L2=importdata([outdir,'LN02_2_line.csv']);
L4=importdata([outdir,'LN02_4_line.csv']);

fds40diamx =strcmp('d32_40-x',L1.colheaders);
fds62diamx =strcmp('d32_60-x',L1.colheaders);
fds40diam  =strcmp('d32_40',L1.colheaders);
fds62diam  =strcmp('d32_60',L1.colheaders);
fds40velox =strcmp('w00_40-x',L1.colheaders);
fds62velox =strcmp('w00_60-x',L1.colheaders);
fds40velo  =strcmp('w00_40',L1.colheaders);
fds62velo  =strcmp('w00_60',L1.colheaders);
fds40fluxx =strcmp('f_40-x',L1.colheaders);
fds62fluxx =strcmp('f_60-x',L1.colheaders);
fds40flux  =strcmp('f_40',L1.colheaders);
fds62flux  =strcmp('f_60',L1.colheaders);

exp40x     =strcmp('Position40',edata.colheaders);
exp62x     =strcmp('Position62',edata.colheaders);
exp40diam  =strcmp('Diam40',edata.colheaders);
exp62diam  =strcmp('Diam62',edata.colheaders);
exp40velo  =strcmp('Velo40',edata.colheaders);
exp62velo  =strcmp('Velo62',edata.colheaders);
exp40flux  =strcmp('Flux40',edata.colheaders);
exp62flux  =strcmp('Flux62',edata.colheaders);


figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X Plot_Y Plot_Width Plot_Height])

plot(edata.data(:,exp40x),-edata.data(:,exp40velo),'ko');
hold on
plot(L1.data(:,fds40velox),-L1.data(:,fds40velo),'k-');
plot(L2.data(:,fds40velox),-L2.data(:,fds40velo),'b-');
plot(L4.data(:,fds40velox),-L4.data(:,fds40velo),'r-');
axis([ -0.3 0.3 0 7]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Velocity (m/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Northeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_velo_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_velo_40');

figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])

plot(edata.data(:,exp40x),edata.data(:,exp40diam)*1e6,'ko');
hold on
plot(L1.data(:,fds40diamx),L1.data(:,fds40diam)*1e6,'k-');
plot(L2.data(:,fds40diamx),L2.data(:,fds40diam)*1e6,'b-');
plot(L4.data(:,fds40diamx),L4.data(:,fds40diam)*1e6,'r-');
axis([ -0.3 0.3 0 160]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Diameter d_{32} ({\mu}m)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Southeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_diam_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_diam_40');

figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])

plot(edata.data(:,exp40x),edata.data(:,exp40flux),'ko');
hold on
plot(L1.data(:,fds40fluxx),-L1.data(:,fds40flux),'k-');
plot(L2.data(:,fds40fluxx),-L2.data(:,fds40flux),'b-');
plot(L4.data(:,fds40fluxx),-L4.data(:,fds40flux),'r-');
axis([ -0.3 0.3 0 0.4]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Droplet Flux (kg/m^2/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Northeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_diam_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_flux_40');


figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])

plot(edata.data(:,exp62x),-edata.data(:,exp62velo),'ko');
hold on
plot(L1.data(:,fds62velox),-L1.data(:,fds62velo),'k-');
plot(L2.data(:,fds62velox),-L2.data(:,fds62velo),'b-');
plot(L4.data(:,fds62velox),-L4.data(:,fds62velo),'r-');
axis([ -0.3 0.3 0 7]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Velocity (m/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Northeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_velo_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_velo_62');

figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])

plot(edata.data(:,exp62x),edata.data(:,exp62diam)*1e6,'ko');
hold on
plot(L1.data(:,fds62diamx),L1.data(:,fds62diam)*1e6,'k-');
plot(L2.data(:,fds62diamx),L2.data(:,fds62diam)*1e6,'b-');
plot(L4.data(:,fds62diamx),L4.data(:,fds62diam)*1e6,'r-');
axis([ -0.3 0.3 0 160]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Diameter d_{32} ({\mu}m)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Southeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_diam_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_diam_62');

figure
set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])

plot(edata.data(:,exp62x),edata.data(:,exp62flux),'ko');
hold on
plot(L1.data(:,fds62fluxx),-L1.data(:,fds62flux),'k-');
plot(L2.data(:,fds62fluxx),-L2.data(:,fds62flux),'b-');
plot(L4.data(:,fds62fluxx),-L4.data(:,fds62flux),'r-');
axis([ -0.3 0.3 0 0.4]);
set(gca,'XTick',[-.3:.1:.3]);
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Radial position (m)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Droplet Flux (kg/m^2/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1 cm','FDS 2 cm','FDS 4 cm'}, 'Location', 'Northeast');
set(h,'Interpreter', Font_Interpreter)
set(h,'FontSize', Key_Font_Size)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'Units', Paper_Units);
set(gcf, 'PaperUnits',Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'Position', [0 0 Paper_Width Paper_Height]);

git_file = [outdir, 'LN02_4_git.txt'];
addverstr(gca,git_file,'linear')

display('Printing plot LN02_diam_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_flux_62');


close all
