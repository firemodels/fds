% McDermott
% 7-23-2009
% beyler_hood.m

close all
clear all

n = 4;

% configuration
alpha{1} = 'CO2';  alt_name{1} = 'CO$_2$'; print_name{1} = 'CO2'; xmax(1) = 0.1;
alpha{2} = 'CO';   alt_name{2} = 'CO';     print_name{2} = 'CO';  xmax(2) = 0.035;
alpha{3} = 'O2';   alt_name{3} = 'O$_2$';  print_name{3} = 'O2';  xmax(3) = 0.16;
alpha{4} = 'C3H8'; alt_name{4} = 'UHC';    print_name{4} = 'UHC'; xmax(4) = 0.07;

addpath('../../Validation/Beyler_Hood/Experimental_Data')
addpath('../../Validation/Beyler_Hood/FDS_Output_Files')

% load experimental data and FDS prediction
[exp_header,exp_data] = dvcread('Beyler_Hood_Test.csv',1);
[fds_header,fds_data] = dvcread('Beyler_Hood_FDS.csv',1);

for j = 1:n

    % plot the results
    i1 = strmatch(['Exp ',alpha{j},' (5 cm)'],exp_header);
    j1 = strmatch(['FDS ',alpha{j},' (5 cm)'],fds_header);
    K(1) = plot(exp_data(:,i1),fds_data(:,j1),'bd','MarkerFaceColor','b'); hold on
    
    i2 = strmatch(['Exp ',alpha{j},' (-10 cm)'],exp_header);
    j2 = strmatch(['FDS ',alpha{j},' (-10 cm)'],fds_header);
    K(2) = plot(exp_data(:,i2),fds_data(:,j2),'rsq','MarkerFaceColor','r');
    
    i3 = strmatch(['Exp ',alpha{j},' (0 cm)'],exp_header);
    j3 = strmatch(['FDS ',alpha{j},' (0 cm)'],fds_header);
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
    ymin = 0;
    ymax = xmax(j);
    
    plot([xmin xmax(j)],[ymin ymax],'k-.')
    axis([xmin xmax(j) ymin ymax])
    xlabel(['Measured ',alt_name{j},' Volume Fraction (\%)'],'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel(['Predicted ',alt_name{j},' Volume Fraction (\%)'],'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    
    h = legend(K,'5 cm','-10 cm','0 cm','Location','NorthWest');
    set(h,'Interpreter',Font_Interpreter)
    
    % add SVN if file is available
    
    SVN_Filename = 'Beyler_Hood_propane-1353-19-0_svn.txt';
    if exist(SVN_Filename,'file')
        SVN = importdata(SVN_Filename);
        x_lim = get(gca,'XLim');
        y_lim = get(gca,'YLim');
        X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
        Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
        text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
            'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
    end
    
    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Height Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Height Paper_Height]);
    print(gcf,'-dpdf',['../../Manuals/FDS_Validation_Guide/FIGURES/Beyler_Hood/Beyler_',print_name{j}])
    
    close
    
end
