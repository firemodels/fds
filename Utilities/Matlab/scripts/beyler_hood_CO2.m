% McDermott
% 7-22-2009
% beyler_hood_CO2.m

close all
clearvars -except repository

%repository = '../../../FDS-SMV';

addpath([repository,'/Validation/Beyler_Hood/Experimental_Data'])
addpath([repository,'/Validation/Beyler_Hood/FDS_Output_Files'])

% load experimental data and FDS prediction
[exp_data,exp_header] = xlsread('Beyler_Hood_Test.csv');
[fds_data,fds_header] = xlsread('Beyler_Hood_FDS.csv');

% plot the results
i1 = strmatch('Exp CO2 (5 cm)',exp_header);
j1 = strmatch('FDS CO2 (5 cm)',fds_header);
K(1) = plot(exp_data(:,i1),fds_data(:,j1),'bd','MarkerFaceColor','b'); hold on

i2 = strmatch('Exp CO2 (-10 cm)',exp_header);
j2 = strmatch('FDS CO2 (-10 cm)',fds_header);
K(2) = plot(exp_data(:,i2),fds_data(:,j2),'rsq','MarkerFaceColor','r');

i3 = strmatch('Exp CO2 (0 cm)',exp_header);
j3 = strmatch('FDS CO2 (0 cm)',fds_header);
K(3) = plot(exp_data(:,i3),fds_data(:,j3),'^g','MarkerFaceColor','g');

% format the plot
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2; 
set(gca,'Position',[Plot_X,Plot_Y,Plot_Height,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

xmin = 0;
xmax = 0.1;
ymin = 0;
ymax = 0.1;
plot([xmin xmax],[ymin ymax],'k-.')
axis([xmin xmax ymin ymax])
xlabel('Measured CO$_2$ Volume Fraction (\%)','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('Predicted CO$_2$ Volume Fraction (\%)','Interpreter','LaTeX','FontSize',Label_Font_Size)

h = legend(K,'5 cm','-10 cm','0 cm','Location','NorthWest');
set(h,'Interpreter','LaTeX')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Height Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Height Paper_Height]);
print(gcf,'-dpdf',[repository,'/Manuals/FDS_5_Validation_Guide/FIGURES/Beyler_Hood/Beyler_CO2'])

