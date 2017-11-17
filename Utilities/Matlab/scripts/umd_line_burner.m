% McDermott
% 4-15-15
% umd_line_burner

close all
clear all

plot_style
Marker_Size = 10;

% % compute N2 ramp

% format long

% MF_AIR = 0.2446; % kg/m2/s
% X_O2 = linspace(2.077972E-01,0.10,61); % White et al. Fig. 5
% X_O2_AIR = 2.077972E-01;
% W_AIR = 28.76431;
% W_N2 = 28.01340;

% % X_O2 = (MF_AIR/W_AIR)*X_O2_AIR / ( MF_AIR/W_AIR + MF_N2/W_N2  )  ---> solve for MF_N2

% MF_N2 = W_N2*( (MF_AIR/W_AIR)*X_O2_AIR./X_O2 - (MF_AIR/W_AIR) );

% F = MF_N2/(MF_N2(end));

% ramp_time = @(x) (x-X_O2(1))*60./(0.1-X_O2(1));
% [round(ramp_time(X_O2))',F']

% return

% compute DEVC positions

% format long

dy1 = 1.25e-2;
yc1 = -(.25-dy1/2):dy1:(.25-dy1/2);

dy2 = .625e-2; % m
yc2 = -(.25-dy2/2):dy2:(.25-dy2/2);

dy3 = .3125e-2; % m
yc3 = -(.25-dy3/2):dy3:(.25-dy3/2);

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.125, QUANTITY=''TEMPERATURE'', TIME_AVERAGED=.FALSE. /'])
% end

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.125, QUANTITY=''THERMOCOUPLE'', PROP_ID=''TC'', TIME_AVERAGED=.FALSE. /'])
% end

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.125, QUANTITY=''VOLUME FRACTION'', SPEC_ID=''OXYGEN'', TIME_AVERAGED=.FALSE. /'])
% end

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.250, QUANTITY=''TEMPERATURE'', TIME_AVERAGED=.FALSE. /'])
% end

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.250, QUANTITY=''THERMOCOUPLE'', PROP_ID=''TC'', TIME_AVERAGED=.FALSE. /'])
% end

% for i=1:length(yc3)
%     display(['&DEVC XYZ=0,',num2str(yc3(i)),',0.250, QUANTITY=''VOLUME FRACTION'', SPEC_ID=''OXYGEN'', TIME_AVERAGED=.FALSE. /'])
% end

% return

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

F1 = importdata([outdir,'methane_dx_1p25cm_line.csv'],',',2);
F2 = importdata([outdir,'methane_dx_p625cm_line.csv'],',',2);
F3 = importdata([outdir,'methane_dx_p3125cm_line.csv'],',',2);

% % oxygen level

% % J1 = importdata([outdir,'methane_dx_1p25cm_devc.csv'],',',2);

% t = J1.data(:,1);
% X_O2 = J1.data(:,2);

% figure
% plot_style

% H(1)=plot(t,X_O2,'k-');

% xlabel('Time (s)')
% ylabel('O2 (vol frac)')

% print(gcf,'-dpdf',[pltdir,'umd_line_burner_O2'])

% return

% % check heat release rate

% HOC = 49674; % kJ/kg methane
% mf = 0.04;   % kg/m2/s methane
% A  = 0.05*0.5; % m2
% G1 = importdata([outdir,'methane_dx_1p25cm_hrr.csv'],',',2);
% t = G1.data(:,1);
% HRR = G1.data(:,strcmp(G1.colheaders,'HRR'));
% MLR = G1.data(:,strcmp(G1.colheaders,'MLR_FUEL'));

% plot(t,ones(length(t),1)*mf*A*HOC,'k--'); hold on
% plot(t,HRR,'r-');
% plot(t,MLR*HOC,'b-');

% xlabel('Time (s)')
% ylabel('HRR (kW)')

% return

% % mass flow for propane case
% HOC = 46334.6246; % kJ/kg propane
% HRRPUA = 2000;
% A = 0.05*.5;
% HRR = HRRPUA*A % should be 50 kW
% mf = 2000/HOC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean temperature at z=0.125 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M1 = importdata([expdir,'TC_Data.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'x_125')));
T1 = M1.data(:,find(strcmp(M1.colheaders,'TC_125')));

H(1)=plot(y1,T1,'ksq','MarkerSize',Marker_Size); hold on

% z = 0.125 m

yc1 = F1.data(:,find(strcmp(F1.colheaders,'TC_125-y')));
TC1 = F1.data(:,find(strcmp(F1.colheaders,'TC_125')));

yc2 = F2.data(:,find(strcmp(F2.colheaders,'TC_125-y')));
TC2 = F2.data(:,find(strcmp(F2.colheaders,'TC_125')));

yc3 = F3.data(:,find(strcmp(F3.colheaders,'TC_125-y')));
TC3 = F3.data(:,find(strcmp(F3.colheaders,'TC_125')));

H(2) = plot(yc1,TC1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,TC2,'m--','LineWidth',Line_Width); % dx = 0.625 cm
H(4) = plot(yc3,TC3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

% % write data to a table file
% Tbl = table([yc1',TC1']);
% writetable(Tbl,[pltdir,'methane_O2_p18_TC_z_p125_table.csv'])

xmin = -0.25;
xmax = 0.25;
ymin = 0;
ymax = 1200;
xt = xmin + .03*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,'UMD Line Burner, CH4','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xt = xmin + .03*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,'18 % O2, {\it z} = 0.125 m','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

axis([xmin xmax ymin ymax])
xlabel('Position (m)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Thermocouple Temperature ( \circC )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
lh = legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','FDS 0.3125 cm');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
git_file = [outdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_TC_z_p125'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean temperature at z=0.250 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M1 = importdata([expdir,'TC_Data.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'x_250')));
T1 = M1.data(:,find(strcmp(M1.colheaders,'TC_250')));

H(1)=plot(y1,T1,'ksq','MarkerSize',Marker_Size); hold on

% z = 0.250 m

yc1 = F1.data(:,find(strcmp(F1.colheaders,'TC_250-y')));
TC1 = F1.data(:,find(strcmp(F1.colheaders,'TC_250')));

yc2 = F2.data(:,find(strcmp(F2.colheaders,'TC_250-y')));
TC2 = F2.data(:,find(strcmp(F2.colheaders,'TC_250')));

yc3 = F3.data(:,find(strcmp(F3.colheaders,'TC_250-y')));
TC3 = F3.data(:,find(strcmp(F3.colheaders,'TC_250')));

H(2) = plot(yc1,TC1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,TC2,'m--','LineWidth',Line_Width); % dx = 0.625 cm
H(4) = plot(yc3,TC3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

% % write data to a table file
% Tbl = table([yc1',TC1']);
% writetable(Tbl,[pltdir,'methane_O2_p18_TC_z_p250_table.csv'])

xmin = -0.25;
xmax = 0.25;
ymin = 0;
ymax = 1200;
xt = xmin + .03*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,'UMD Line Burner, CH4','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xt = xmin + .03*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,'18 % O2, {\it z} = 0.250 m','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

axis([xmin xmax ymin ymax])
xlabel('Position (m)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Thermocouple Temperature ( \circC )','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
lh=legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','FDS 0.3125 cm');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
git_file = [outdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_TC_z_p250'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O2 concentration at z=0.125 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M1 = importdata([expdir,'O2_Data.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'x_125')));
O21 = M1.data(:,find(strcmp(M1.colheaders,'XO2_125')));

H(1)=plot(y1,O21,'ko','MarkerSize',Marker_Size); hold on

% z = 0.125 m

yc1 = F1.data(:,find(strcmp(F1.colheaders,'XO2_125-y')));
O2_1 = F1.data(:,find(strcmp(F1.colheaders,'XO2_125')));

yc2 = F2.data(:,find(strcmp(F2.colheaders,'XO2_125-y')));
O2_2 = F2.data(:,find(strcmp(F2.colheaders,'XO2_125')));

yc3 = F3.data(:,find(strcmp(F3.colheaders,'XO2_125-y')));
O2_3 = F3.data(:,find(strcmp(F3.colheaders,'XO2_125')));

H(2) = plot(yc1,O2_1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,O2_2,'m--','LineWidth',Line_Width); % dx = 0.625 cm
H(4) = plot(yc3,O2_3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

% % write data to a table file
% Tbl = table([yc1',O2_1']);
% writetable(Tbl,[pltdir,'methane_O2_p18_O2_z_p125_table.csv'])

xmin = -0.25;
xmax = 0.25;
ymin = 0.05;
ymax = 0.25;
xt = xmin + .03*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,'UMD Line Burner, CH4','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xt = xmin + .03*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,'18 % O2, {\it z} = 0.125 m','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

axis([xmin xmax ymin ymax])
xlabel('Position (m)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('O2 (vol frac)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
lh=legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','FDS 0.3125 cm','Location','Southwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
git_file = [outdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_O2_z_p125'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O2 concentration at z=0.250 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M1 = importdata([expdir,'O2_Data.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'x_250')));
O21 = M1.data(:,find(strcmp(M1.colheaders,'XO2_250')));

H(1)=plot(y1,O21,'ko','MarkerSize',Marker_Size); hold on

% z = 0.250 m

yc1 = F1.data(:,find(strcmp(F1.colheaders,'XO2_250-y')));
O2_1 = F1.data(:,find(strcmp(F1.colheaders,'XO2_250')));

yc2 = F2.data(:,find(strcmp(F2.colheaders,'XO2_250-y')));
O2_2 = F2.data(:,find(strcmp(F2.colheaders,'XO2_250')));

yc3 = F3.data(:,find(strcmp(F3.colheaders,'XO2_250-y')));
O2_3 = F3.data(:,find(strcmp(F3.colheaders,'XO2_250')));

H(2) = plot(yc1,O2_1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,O2_2,'m--','LineWidth',Line_Width);  % dx = 0.625 cm
H(4) = plot(yc3,O2_3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

% % write data to a table file
% Tbl = table([yc1',O2_1']);
% writetable(Tbl,[pltdir,'methane_O2_p18_O2_z_p250_table.csv'])

xmin = -0.25;
xmax = 0.25;
ymin = 0.05;
ymax = 0.25;
xt = xmin + .03*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,'UMD Line Burner, CH4','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xt = xmin + .03*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,'18 % O2, {\it z} = 0.250 m','FontName',Font_Name,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

axis([xmin xmax ymin ymax])
xlabel('Position (m)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('O2 (vol frac)','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
lh=legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','FDS 0.3125 cm','Location','Southwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
git_file = [outdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_O2_z_p250'])
