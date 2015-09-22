% McDermott
% 4-15-15
% umd_line_burner

close all
clear all

% % compute N2 ramp

% format long

% MF_AIR = 0.2446; % kg/m2/s
% X_O2 = [0.210, 0.181, 0.168, 0.158, 0.151]; % White et al. Fig. 5
% X_O2_AIR = 2.077972E-01;
% W_AIR = 28.76431;
% W_N2 = 28.01340;

% % X_O2 = (MF_AIR/W_AIR)*X_O2_AIR / ( MF_AIR/W_AIR + MF_N2/W_N2  )  ---> solve for MF_N2

% MF_N2 = W_N2*( (MF_AIR/W_AIR)*X_O2_AIR./X_O2 - (MF_AIR/W_AIR) )

% F = MF_N2/(MF_N2(end))

% return

% % compute DEVC positions

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

expdir = '../../Validation/UMD_Line_Burner/Experimental_Data/';
fdsdir = '../../Validation/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

F1 = importdata([fdsdir,'methane_dx_1p25cm_devc.csv'],',',2);
F2 = importdata([fdsdir,'methane_dx_p625cm_devc.csv'],',',2);
%F3 = importdata([fdsdir,'methane_dx_p3125cm_devc.csv'],',',2);

% % oxygen level

% % J1 = importdata([fdsdir,'methane_extinction_dx_1p25cm_devc.csv'],',',2);

% t = J1.data(:,1);
% X_O2 = J1.data(:,2);

% figure
% plot_style

% H(1)=plot(t,X_O2,'k-');

% xlabel('Time (s)')
% ylabel('O2 (vol frac)')

% print(gcf,'-dpdf',[pltdir,'umd_gas_burner_O2'])

% return

% % check heat release rate

% HOC = 49674; % kJ/kg methane
% mf = 0.04;   % kg/m2/s methane
% A  = 0.05*0.5; % m2
% G1 = importdata([fdsdir,'methane_dx_1p25cm_hrr.csv'],',',2);
% t = G1.data(:,1);
% HRR = G1.data(:,strcmp(G1.colheaders,'HRR'));
% MLR = G1.data(:,strcmp(G1.colheaders,'MLR_FUEL'));

% plot(t,ones(length(t),1)*mf*A*HOC,'k--'); hold on
% plot(t,HRR,'r-');
% plot(t,MLR*HOC,'b-');

% xlabel('Time (s)')
% ylabel('HRR (kW)')

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean temperature at z=0.125 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot_style
Marker_Size = 10;

M1 = importdata([expdir,'Exp_O2_p18_T_z_p125m.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'y (m)')));
T1 = M1.data(:,find(strcmp(M1.colheaders,'T (C)')));

H(1)=plot(y1,T1,'ksq','MarkerSize',Marker_Size); hold on

% find 18 % O2 time block

t1 = F1.data(:,1); t1range = find(t1>=30 & t1<=40);
t2 = F2.data(:,1); t2range = find(t2>=30 & t2<=40);
%t3 = F3.data(:,1); t3range = find(t3>=30 & t3<=40);

C1 = 3;
ND1 = 40;
TG1_p125 = C1:C1+ND1-1;
TC1_p125 = C1+ND1:C1+2*ND1-1;

C2 = 3;
ND2 = 80;
TG2_p125 = C2:C2+ND2-1;
TC2_p125 = C2+ND2:C2+2*ND2-1;

% C3 = 3;
% ND3 = 160;
% TG3_p125 = C3:C3+ND3-1;
% TC3_p125 = C3+ND3:C3+2*ND3-1;

% z = 0.125 m

TG1 = mean(F1.data(t1range,TG1_p125),1);
TC1 = mean(F1.data(t1range,TC1_p125),1);

TG2 = mean(F2.data(t2range,TG2_p125),1);
TC2 = mean(F2.data(t2range,TC2_p125),1);

% TG3 = mean(F3.data(t3range,TG3_p125),1);
% TC3 = mean(F3.data(t3range,TC3_p125),1);

H(2) = plot(yc1,TC1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,TC2,'m--','LineWidth',Line_Width);  % dx = 0.625 cm
%H(4) = plot(yc3,TC3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

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
xlabel('Position (m)')
ylabel('Thermocouple Temperature ( \circC )')
lh = legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm'); %,'FDS 0.3125 cm');
set(lh,'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
git_file = [fdsdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_TC_z_p125'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean temperature at z=0.250 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot_style
Marker_Size = 10;

M1 = importdata([expdir,'Exp_O2_p18_T_z_p250m.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'y (m)')));
T1 = M1.data(:,find(strcmp(M1.colheaders,'T (C)')));

H(1)=plot(y1,T1,'ksq','MarkerSize',Marker_Size); hold on

% find 18 % O2 time block

t1 = F1.data(:,1); t1range = find(t1>=30 & t1<=40);
t2 = F2.data(:,1); t2range = find(t2>=30 & t2<=40);
%t3 = F3.data(:,1); t3range = find(t3>=30 & t3<=40);

C1 = 3;
ND1 = 40;
TG1_p250 = C1+3*ND1:C1+4*ND1-1;
TC1_p250 = C1+4*ND1:C1+5*ND1-1;

C2 = 3;
ND2 = 80;
TG2_p250 = C2+3*ND2:C2+4*ND2-1;
TC2_p250 = C2+4*ND2:C2+5*ND2-1;

% C3 = 3;
% ND3 = 160;
% TG3_p250 = C3+3*ND3:C3+4*ND3-1;
% TC3_p250 = C3+4*ND3:C3+5*ND3-1;

TG1 = mean(F1.data(t1range,TG1_p250),1);
TC1 = mean(F1.data(t1range,TC1_p250),1);

TG2 = mean(F2.data(t2range,TG2_p250),1);
TC2 = mean(F2.data(t2range,TC2_p250),1);

% TG3 = mean(F3.data(t3range,TG3_p250),1);
% TC3 = mean(F3.data(t3range,TC3_p250),1);

H(2) = plot(yc1,TC1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,TC2,'m--','LineWidth',Line_Width);  % dx = 0.625 cm
%H(4) = plot(yc3,TC3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

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
xlabel('Position (m)')
ylabel('Thermocouple Temperature ( \circC )')
legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm') %,'FDS 0.3125 cm')

% add Git revision if file is available
git_file = [fdsdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_TC_z_p250'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O2 concentration at z=0.125 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot_style
Marker_Size = 10;

M1 = importdata([expdir,'Exp_O2_p18_O2_z_p125m.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'y (m)')));
O21 = M1.data(:,find(strcmp(M1.colheaders,'O2 (vol frac)')));

H(1)=plot(y1,O21,'ko','MarkerSize',Marker_Size); hold on

% find 18 % O2 time block

t1 = F1.data(:,1); t1range = find(t1>=30 & t1<=40);
t2 = F2.data(:,1); t2range = find(t2>=30 & t2<=40);
%t3 = F3.data(:,1); t3range = find(t3>=30 & t3<=40);

C1 = 3;
ND1 = 40;
O21_p125 = C1+2*ND1:C1+3*ND1-1;

C2 = 3;
ND2 = 80;
O22_p125 = C2+2*ND2:C2+3*ND2-1;

% C3 = 3;
% ND3 = 160;
% O23_p125 = C3+2*ND3:C3+3*ND3-1;

% z = 0.125 m

O2_1 = mean(F1.data(t1range,O21_p125),1);
O2_2 = mean(F2.data(t2range,O22_p125),1);
%O2_3 = mean(F3.data(t3range,O23_p125),1);

H(2) = plot(yc1,O2_1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,O2_2,'m--','LineWidth',Line_Width);  % dx = 0.625 cm
%H(4) = plot(yc3,O2_3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

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
xlabel('Position (m)')
ylabel('O2 (vol frac)')
legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','Location','Southwest') %,'FDS 0.3125 cm','Location','Southwest')

% add Git revision if file is available
git_file = [fdsdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_O2_z_p125'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O2 concentration at z=0.250 m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot_style
Marker_Size = 10;

M1 = importdata([expdir,'Exp_O2_p18_O2_z_p250m.csv'],',',1);

y1 = M1.data(:,find(strcmp(M1.colheaders,'y (m)')));
O21 = M1.data(:,find(strcmp(M1.colheaders,'O2 (vol frac)')));

H(1)=plot(y1,O21,'ko','MarkerSize',Marker_Size); hold on

% find 18 % O2 time block

t1 = F1.data(:,1); t1range = find(t1>=30 & t1<=40);
t2 = F2.data(:,1); t2range = find(t2>=30 & t2<=40);
%t3 = F3.data(:,1); t3range = find(t3>=30 & t3<=40);

C1 = 3;
ND1 = 40;
O21_p250 = C1+5*ND1:C1+6*ND1-1;

C2 = 3;
ND2 = 80;
O22_p250 = C2+5*ND2:C2+6*ND2-1;

% C3 = 3;
% ND3 = 160;
% O23_p250 = C3+5*ND3:C3+6*ND3-1;

% z = 0.125 m

O2_1 = mean(F1.data(t1range,O21_p250),1);
O2_2 = mean(F2.data(t2range,O22_p250),1);
%O2_3 = mean(F3.data(t3range,O23_p250),1);

H(2) = plot(yc1,O2_1,'r-.','LineWidth',Line_Width); % dx = 1.25 cm
H(3) = plot(yc2,O2_2,'m--','LineWidth',Line_Width);  % dx = 0.625 cm
%H(4) = plot(yc3,O2_3,'b-','LineWidth',Line_Width);  % dx = 0.3125 cm

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
xlabel('Position (m)')
ylabel('O2 (vol frac)')
legend(H,'Exp','FDS 1.25 cm','FDS 0.625 cm','Location','Southwest') %,'FDS 0.3125 cm','Location','Southwest')

% add Git revision if file is available
git_file = [fdsdir,'methane_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[pltdir,'methane_O2_p18_O2_z_p250'])