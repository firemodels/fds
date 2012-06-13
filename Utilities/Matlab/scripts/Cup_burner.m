% Hostikka
% 3-30-2010
% cup_burner.m

close all
clear all

addpath('../../Validation/Cup_Burner/Experimental_Data')
addpath('../../Validation/Cup_Burner/FDS_Output_Files')

% load experimental data and FDS prediction
[exp_data] = csvread('Cup_burner_data.csv',2);

HRR_Limit = 1E-10;

N_Fuels = 2;
Fuel{1} = 'CH4';
Fuel{2} = 'C7H16';

N_Agents = 4;
Agent{1} = 'He';
Agent{2} = 'Ar';
Agent{3} = 'N2';
Agent{4} = 'CO2';

X_leg_pos = [0.55 0.3 0.2 0.2];
Y_leg_pos = [0.55 0.3 0.2 0.2];

% Color per fuel
color{1} = 'b';
color{2} = 'r';

% Marker per Agens
marker{1} = 'o';
marker{2} = 'x';
marker{3} = '+';
marker{4} = '*';
marker{5} = 's';
marker{6} = 'd';
marker{7} = '^';
marker{8} = '<';
marker{9} = '>';

MarkerSize = 7;
LineWidth = 1.5;

% Collect data

for f = 1:N_Fuels
for s = 1:N_Agents

   FDS_File = ['Cup_' Fuel{f} '_' Agent{s} '_devc.csv'];
   [fds_data] = csvread(FDS_File,2);
   n_fds = size(fds_data,1);
   i_first = find(fds_data(1:n_fds,4)>HRR_Limit,1,'first');
   i_last = find(fds_data(i_first:n_fds,4)<HRR_Limit,1,'first')+(i_first-1);
   if (isempty(i_last))
      % not yet extinguished
      i_last = 1;
   end
   FDS_X(f,s) = fds_data(i_last,3);
   FDS_Y(f,s) = fds_data(i_last,2);

   Exp_X(f,s) = exp_data(1,(f-1)*N_Agents+s);
   Exp_Y(f,s) = exp_data(2,(f-1)*N_Agents+s);   
   
end
end

Xmax = max(max(FDS_X));
Xmax = max(max(max(Exp_X)),Xmax);
Xmax = ceil(Xmax*10)/10;
Ymax = max(max(FDS_Y));
Ymax = max(max(max(Exp_Y)),Ymax);
Ymax = ceil(Ymax*10)/10;

% plot X

hf(1)=figure(1);
n = 0;
for f = 1:N_Fuels
for s = 1:N_Agents
   n = n + 1;
   hX(n) = plot(Exp_X(f,s),FDS_X(f,s));
   XLegendStr{n} = [Fuel{f} ' ' Agent{s}];
   set(hX(n),'Marker',marker{s},...
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
xlabel('Measured MEC (volume fraction)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel('Predicted MEC (volume fraction)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
%legend(hX,XLegendStr,'Location',X_leg_pos)
legend(hX,XLegendStr,'Location','NorthWest')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Validation_Guide/FIGURES/Cup_Burner/Cup_Burner_volfrac');


% plot Y

hf(2)=figure(2);
n = 0;
for f = 1:N_Fuels
for s = 1:N_Agents
   n = n + 1;
   hY(n) = plot(Exp_Y(f,s),FDS_Y(f,s));
   YLegendStr{n} = [Fuel{f} ' ' Agent{s}];
   set(hY(n),'Marker',marker{s},...
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
xmax = Ymax;
ymax = xmax;
plot([xmin xmax],[ymin ymax],'k-.')
axis([xmin xmax ymin ymax])

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2;
set(gca,'Position',[Plot_X,Plot_Y,Plot_Height,Plot_Height])
set(hf(2),'DefaultLineLineWidth',Line_Width)
xlabel('Measured MEC (mass fraction)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel('Predicted MEC (mass fraction)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
%legend(hY,YLegendStr,'Location',Y_leg_pos)
legend(hY,YLegendStr,'Location','NorthWest')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Validation_Guide/FIGURES/Cup_Burner/Cup_Burner_massfrac');

%close all