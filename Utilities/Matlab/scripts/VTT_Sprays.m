% Topi Sikanen

plot_style

repository = '../../Validation/VTT_Sprays/';

files = dir(strcat(repository,'*_line.csv'));

chids=strrep({files.name},'FDS_Output_Files/*_line.csv','');
cases={ 'LN02_1' ,'LN02_2','LN02_4'};
stop = ~all(cellfun(@(cname) any(strcmp(cname,cases)),chids));

if(stop) 
    display(['Error: File missing in ' repository, '. Skipping case.']);
    return
end
   
edata=importdata([repository,'Experimental_Data/LN02.csv']);
L1=importdata([repository,'FDS_Output_Files/LN02_1_line.csv']);
L2=importdata([repository,'FDS_Output_Files/LN02_2_line.csv']);
L4=importdata([repository,'FDS_Output_Files/LN02_4_line.csv']);



fds40diamx =strcmp('d10_40-x',L1.colheaders);
fds62diamx =strcmp('d10_60-x',L1.colheaders);
fds40diam  =strcmp('d10_40',L1.colheaders);
fds62diam  =strcmp('d10_60',L1.colheaders);
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


figure(1)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp40x),-edata.data(:,exp40velo),'ko');
hold on
plot(L1.data(:,fds40velox),-L1.data(:,fds40velo),'k-');
plot(L2.data(:,fds40velox),-L2.data(:,fds40velo),'b-');
plot(L4.data(:,fds40velox),-L4.data(:,fds40velo),'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Velocity (m/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Northeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_velo_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_velo_40');

figure(2)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp40x),edata.data(:,exp40diam)*1e6,'ko');
hold on
plot(L1.data(:,fds40diamx),L1.data(:,fds40diam)*1e6,'k-');
plot(L2.data(:,fds40diamx),L2.data(:,fds40diam)*1e6,'b-');
plot(L4.data(:,fds40diamx),L4.data(:,fds40diam)*1e6,'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Diameter ({\mu}m)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Southeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_diam_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_diam_40');

figure(3)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp40x),edata.data(:,exp40flux),'ko');
hold on
plot(L1.data(:,fds40fluxx),-L1.data(:,fds40flux),'k-');
plot(L2.data(:,fds40fluxx),-L2.data(:,fds40flux),'b-');
plot(L4.data(:,fds40fluxx),-L4.data(:,fds40flux),'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Droplet Flux (kg/m^2s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Northeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_diam_40.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_flux_40');


figure(4)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp62x),-edata.data(:,exp62velo),'ko');
hold on
plot(L1.data(:,fds62velox),-L1.data(:,fds62velo),'k-');
plot(L2.data(:,fds62velox),-L2.data(:,fds62velo),'b-');
plot(L4.data(:,fds62velox),-L4.data(:,fds62velo),'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Velocity (m/s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Northeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_velo_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_velo_62');

figure(5)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp62x),edata.data(:,exp62diam)*1e6,'ko');
hold on
plot(L1.data(:,fds62diamx),L1.data(:,fds62diam)*1e6,'k-');
plot(L2.data(:,fds62diamx),L2.data(:,fds62diam)*1e6,'b-');
plot(L4.data(:,fds62diamx),L4.data(:,fds62diam)*1e6,'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Diameter ({\mu}m)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Northeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_diam_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_diam_62');

figure(6)

set(gca, 'Units', Plot_Units)
set(gca, 'Position', [Plot_X, Plot_Y, Plot_Width, Plot_Height])
set(gcf, 'DefaultLineLineWidth', Line_Width)

plot(edata.data(:,exp62x),edata.data(:,exp62flux),'ko');
hold on
plot(L1.data(:,fds62fluxx),-L1.data(:,fds62flux),'k-');
plot(L2.data(:,fds62fluxx),-L2.data(:,fds62flux),'b-');
plot(L4.data(:,fds62fluxx),-L4.data(:,fds62flux),'r-');

set(gca, 'FontName', Font_Name)
set(gca, 'FontSize', Key_Font_Size)

xlabel('Radial position (cm)', 'Interpreter', Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Droplet Flux (kg/m^2s)','FontSize',Label_Font_Size)
h = legend({'Experiment','FDS 1cm','FDS 2cm','FDS 4cm'}, 'Location', 'Northeast');

set(h,'Interpreter', Font_Interpreter)

set(gcf, 'Visible', Figure_Visibility);
set(gcf, 'PaperUnits', Paper_Units);
set(gcf, 'PaperSize', [Paper_Width Paper_Height]);
set(gcf, 'PaperPosition', [0 0 Paper_Width Paper_Height]);

SVN_Filename = [repository, 'FDS_Output_Files/LN02_4_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

display('Printing plot LN02_diam_62.pdf...')
print(gcf, '-dpdf', '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/VTT_Sprays/LN02_flux_62');


close all