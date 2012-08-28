% Floyd
% 8-28-2012
% Beyler_hood.m

close all
clear all

addpath('../../Validation/Beyler_Hood/Experimental_Data')
addpath('../../Validation/Beyler_Hood/FDS_Output_Files')

% load experimental data and FDS prediction
[exp_data] = csvread('Beyler_Hood_data.csv',2);

Fuel{1} = 'acetone';

TestID(1,1)  = 117;
TestID(1,2)  = 118;
TestID(1,3)  = 119;
TestID(1,4)  = 120;
TestID(1,5)  = 121;
TestID(1,6)  = 122;
TestID(1,7)  = 124;
TestID(1,8)  = 125;
TestID(1,9)  = 126;
TestID(1,10) = 128;
TestID(1,11) = 143;
TestID(1,12) = 143;
TestID(1,13) = 145;

Fuel{2} = 'ethanol';

TestID(2,1)  = 106;
TestID(2,2)  = 107;
TestID(2,3)  = 108;
TestID(2,4)  = 109;
TestID(2,5)  = 110;
TestID(2,6)  = 111;
TestID(2,7)  = 112;
TestID(2,8)  = 113;
TestID(2,9)  = 114;
TestID(2,10) = 115;
TestID(2,11) = 116;

Fuel{3} = 'isopropanol';

TestID(3,1)  = 129;
TestID(3,2)  = 130;
TestID(3,3)  = 131;
TestID(3,4)  = 132;
TestID(3,5)  = 133;
TestID(3,6)  = 135;
TestID(3,7)  = 136;
TestID(3,8)  = 137;
TestID(3,9)  = 138;
TestID(3,10) = 140;
TestID(3,11) = 141;

Fuel{4} = 'methanol';

TestID(3,1)  = 938;
TestID(3,2)  = 939;
TestID(3,3)  = 940;
TestID(3,4)  = 941;
TestID(3,5)  = 942;
TestID(3,6)  = 943;
TestID(3,7)  = 945;
TestID(3,8)  = 946;
TestID(3,9)  = 947;
TestID(3,10) = 948;
TestID(3,11) = 949;
TestID(3,12) = 950;
TestID(3,13) = 951;
TestID(3,14) = 953;
TestID(3,15) = 954;
TestID(3,16) = 955;
TestID(3,17) = 956;
TestID(3,18) = 957;

Fuel{5} = 'propane';

TestID(5,1)  = 220;
TestID(5,2)  = 228;
TestID(5,3)  = 232;
TestID(5,4)  = 239;
TestID(5,5)  = 241;
TestID(5,6)  = 245;
TestID(5,7)  = 253;
TestID(5,8)  = 257;
TestID(5,9)  = 269;
TestID(5,10) = 278;
TestID(5,11) = 282;
TestID(5,12) = 286;
TestID(5,13) = 287;
TestID(5,14) = 291;
TestID(5,15) = 295;
TestID(5,16) = 299;
TestID(5,17) = 303;
TestID(5,18) = 307;
TestID(5,19) = 311;
TestID(5,20) = 314;
TestID(5,21) = 318;
TestID(5,22) = 322;
TestID(5,23) = 326;
TestID(5,24) = 330;
TestID(5,25) = 334;
TestID(5,26) = 338;
TestID(5,27) = 351;
TestID(5,28) = 355;
TestID(5,29) = 359;
TestID(5,30) = 363;
TestID(5,31) = 367;
TestID(5,32) = 371;
TestID(5,33) = 385;
TestID(5,34) = 389;
TestID(5,35) = 393;
TestID(5,36) = 399;
TestID(5,37) = 403;
TestID(5,38) = 407;
TestID(5,39) = 412;
TestID(5,40) = 417;
TestID(5,41) = 421;
TestID(5,42) = 425;
TestID(5,43) = 429;
TestID(5,44) = 433;
TestID(5,45) = 437;
TestID(5,46) = 441;
TestID(5,47) = 445;


Fuel{6} = 'propylene';

TestID(6,1)  = 776;
TestID(6,2)  = 780;
TestID(6,3)  = 784;
TestID(6,4)  = 792;
TestID(6,5)  = 801;
TestID(6,6)  = 805;
TestID(6,7)  = 809;
TestID(6,8)  = 813;
TestID(6,9)  = 838;
TestID(6,10) = 842;
TestID(6,11) = 846;
TestID(6,12) = 850;
TestID(6,13) = 854;
TestID(6,14) = 859;
TestID(6,15) = 863;
TestID(6,16) = 867;
TestID(6,17) = 870;
TestID(6,18) = 874;
TestID(6,19) = 878;
TestID(6,20) = 882;
TestID(6,21) = 886;
TestID(6,22) = 890;
TestID(6,23) = 899;
TestID(6,24) = 903;
TestID(6,25) = 910;
TestID(6,26) = 914;

Fuel{7} = 'toluene';

TestID(7,1)  = 160;
TestID(7,2)  = 161;
TestID(7,3)  = 162;
TestID(7,4)  = 163;
TestID(7,5)  = 164;
TestID(7,6)  = 165;
TestID(7,7)  = 166;
TestID(7,8)  = 167;
TestID(7,9)  = 168;
TestID(7,10) = 170;

Species{1} = 'O_2';
Species{2} = 'CO_2';
Species{3} = 'H_2O';
Species{4} = 'CO';
Species{5} = 'UHC';
Species{6} = 'Soot';

NumPoints(1) = 13;
NumPoints(2) = 11;
NumPoints(3) = 11;
NumPoints(4) = 18;
NumPoints(5) = 47;
NumPoints(6) = 26;
NumPoints(7) = 10;

N_Fuels = 7;
N_Species = 6;

X_leg_pos = [0.55 0.3 0.2 0.2];
Y_leg_pos = [0.55 0.3 0.2 0.2];

% Color per fuel
color{1} = 'k';
color{2} = 'r';
color{3} = 'b';
color{4} = 'g';
color{5} = 'm';
color{6} = 'c';
color{7} = 'y';

% Marker per fuel
marker{1} = 'o';
marker{2} = 's';
marker{3} = 'd';
marker{4} = '^';
marker{5} = 'v';
marker{6} = '>';
marker{7} = '<';

MarkerSize = 7;
LineWidth = 1.5;

% Collect data

for f = 1:N_Fuels
   for s = 1:NumPoints(f)
      FDS_File = ['Beyler_Hood_' Fuel{f} '_' num2str(TestID(f,s)) '_devc.csv']
      [fds_data] = csvread(FDS_File,2);
      n_fds = size(fds_data,1);
      for ns = 1:N_Species
         ExpPlot(f,s,ns) = exp_data(s,(f-1)*N_Fuels+ns);
         FDSPlot(f,s,ns) = mean(fds_data(n_fds-200:n_fds,253+ns));
      end   
   end
end

for ns = 1:N_Species
   hf(ns)=figure(ns);
   n = 0; 
   Xmax = max(max(FDSPlot(:,:,ns)))
   Xmax = max(max(max(ExpPlot(:,:,ns))),Xmax)
   Xmax = ceil(Xmax*10)/10;
   for s = 1:NumPoints(f)
      for f = 1:N_Fuels
         XLegendStr{f} = [Fuel{f}];
         n = n + 1;
         hX(n) = plot(ExpPlot(f,s,ns),FDSPlot(f,s,ns));
         set(hX(n),'Marker',marker{f},...
        'MarkerSize',MarkerSize,...
        'MarkerEdgeColor',color{f},...
        'MarkerFaceColor','none',...
        'LineWidth',LineWidth,...
        'LineStyle','none');
        hold on
      end
   end
   
xmin = 0;
ymin = 0;
xmax = Xmax;
ymax = xmax;
plot([xmin xmax],[ymin ymax],'k-.')
axis([xmin xmax ymin ymax])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2;
set(gca,'Position',[Plot_X,Plot_Y,Plot_Height,Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xtitle = ['Measured $' Species{ns} '$ (volume fraction)']
ytitle = ['Predicted $' Species{ns} '$ (volume fraction)']
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
%legend(hX,XLegendStr,'Location',X_leg_pos)
legend(hX,XLegendStr,'Location','NorthWest')

% add SVN if file is available

svn_file = '../../Validation/Beyler_Hood/FDS_Output_Files/Beyler_Hood_acetone_117_svn.txt';

if exist(svn_file,'file')
    SVN = importdata(svn_file);
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
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/FIGURES/Beyler_Hood/Beyler_Hood_' Species{ns}]
print(gcf,'-dpdf',plotname);
   
clear hX

end


