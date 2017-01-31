% Hostikka
% 3-30-2010
% cup_burner.m

close all
clear all

expdir = '../../../exp/Cup_Burner/';
outdir = '../../../out/Cup_Burner/FDS_Output_Files/';

plot_style

% load experimental data and FDS prediction
[exp_data] = csvread([expdir,'Cup_burner_data.csv'],2);

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

MarkerSize = 10;

% Collect data

for f = 1:N_Fuels
for s = 1:N_Agents

   FDS_File = [outdir 'Cup_' Fuel{f} '_' Agent{s} '_devc.csv'];
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
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

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
      'LineWidth',Line_Width,...
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

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xlabel('Measured MEC (volume fraction)','FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
ylabel('Predicted MEC (volume fraction)','FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
lh=legend(hX,XLegendStr,'Location','NorthWest');
set(lh,'FontSize',Key_Font_Size)

% add VerStr if file is available

git_file = [outdir,'Cup_C7H16_CO2_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Cup_Burner/Cup_Burner_volfrac');

% plot Y

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

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
      'LineWidth',Line_Width,...
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
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xlabel('Measured MEC (mass fraction)','FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
ylabel('Predicted MEC (mass fraction)','FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
lh=legend(hY,YLegendStr,'Location','NorthWest');
set(lh,'FontSize',Key_Font_Size)

% add VerStr if file is available

git_file = [outdir,'Cup_C7H16_CO2_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Cup_Burner/Cup_Burner_massfrac');

