% Hostikka
% 3-25-2010
% BRE_spray.m

close all
clear all

addpath('../../Validation/BRE_Spray/Experimental_Data')
addpath('../../Validation/BRE_Spray/FDS_Input_Files')
addpath('../../Validation/BRE_Spray/FDS_Output_Files')

% load experimental data and FDS prediction

[exp_data] = csvread('BRE_Spray_Test.csv',2);
exp_col_d = 6;
exp_col_w = 7;
exp_col_att = 4;

amax = max(exp_data(:,exp_col_att));

Nozzle{1} = 'A';
Nozzle{2} = 'B';
Nozzle{3} = 'D';

exp_marker{1} = 'bd';
exp_marker{2} = 'gs';
exp_marker{3} = 'ro';
exp_color{1} = 'b';
exp_color{2} = 'g';
exp_color{3} = 'r';

w_leg_pos{1} = 'NorthWest';
w_leg_pos{2} = 'SouthWest';
w_leg_pos{3} = 'NorthWest';
d_leg_pos{1} = 'NorthWest';
d_leg_pos{2} = 'NorthWest';
d_leg_pos{3} = 'SouthWest';

exp_rows(1,:) = [1:8];
exp_rows(2,:) = [9:16];
exp_rows(3,:) = [17:24];

fds_marker{1} = 'bd--';
fds_marker{2} = 'gs--';
fds_marker{3} = 'ro--';
fds_color{1} = 'w';
fds_color{2} = 'w';
fds_color{3} = 'w';

wmin = 0;
wmax = 6;
dmin = 0;
dmax = 800;

for n = 1:3
for p = 1:8

   FDS_File = ['BRE_Spray_' Nozzle{n} '_' int2str(p) '_devc.csv'];
   [fds_data] = csvread(FDS_File,2);
   
   % Collect data
   i1 = find(fds_data(:,1)>=0.4,1,'first');
   i2 = find(fds_data(:,1)<=0.5,1,'last');
   FDS_Flux_0(n,p) = mean(fds_data(i1:i2,3));
   i1 = find(fds_data(:,1)>=5,1,'first');
   i2 = find(fds_data(:,1)>0.5,1,'last');
   FDS_Flux(n,p) = mean(fds_data(i1:i2,3));
   FDS_d32(n,p) = fds_data(i2,4)*1E6;
   dmax = max(dmax,FDS_d32(n,p));
   FDS_w(n,p) = -1*fds_data(i2,5);
   wmax = max(wmax,FDS_w(n,p));
   FDS_Attenuation(n,p) = 100*(FDS_Flux_0(n,p)-FDS_Flux(n,p))/FDS_Flux_0(n,p);
   amax = max(amax,FDS_Attenuation(n,p));

   
  
end

% plot dv50

hf(1)=figure(1);
subplot(1,3,n)
d_ax(n) = gca;
hd(n,1) = plot(1:8,exp_data(exp_rows(n,:),6),exp_marker{n});
set(hd(n,1),'MarkerFaceColor',exp_color{n});
hold on
hd(n,2) = plot(1:8,FDS_d32(n,:),fds_marker{n});
set(hd(n,2),'MarkerFaceColor',fds_color{n});

% plot w

hf(2)=figure(2);
subplot(1,3,n)
hw(n,1) = plot(1:8,exp_data(exp_rows(n,:),7),exp_marker{n});
w_ax(n) = gca;
set(hw(n,1),'MarkerFaceColor',exp_color{n});
hold on
hw(n,2) = plot(1:8,FDS_w(n,:),fds_marker{n});
set(hw(n,2),'MarkerFaceColor',fds_color{n});

% plot attenuation

hf(3)=figure(3);
ha(n,1) = plot(exp_data(exp_rows(n,:),4),FDS_Attenuation(n,:),exp_marker{n});
set(ha(n,1),'MarkerFaceColor',exp_color{n});
hold on

xmin = 0;
ymin = 0;
xmax = amax;
ymax = amax;
plot([xmin xmax],[ymin ymax],'k-.')
axis([xmin xmax ymin ymax])
  
end


% Format the plots

% d-plot
figure(hf(1))
plot_style
set(d_ax(1:3),'Units',Plot_Units)
set(d_ax(1:3),'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2;
set(gcf,'DefaultLineLineWidth',Line_Width)
for n = 1:3
   set(d_ax(n),'Position',[Plot_X+(n-1)*Plot_Height/3,Plot_Y,Plot_Height/3,Plot_Height])
   axis(d_ax(n),[0 9 dmin dmax]);
   xlabel(d_ax(n),'P (bar)','Interpreter','LaTeX','FontSize',Label_Font_Size)
   title(d_ax(n),['Nozzle ' Nozzle{n}],'Interpreter','LaTeX','FontSize',Label_Font_Size)
   hl(n) = legend(hd(n,:),'Exp.','FDS','Location',d_leg_pos{n});
end
set(d_ax(2:3),'YTickLabel',[])
set(d_ax(1),'XTick',[0 2 4 6 8])
set(d_ax(2:3),'XTick',[2 4 6 8])
ylabel(d_ax(1),'Mean diameter ($\mu$m)','Interpreter','LaTeX','FontSize',Label_Font_Size)

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Height Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Height Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_5_Validation_Guide/FIGURES/BRE_Spray/BRE_Spray_Diameter');


% w-plot

figure(hf(2))
plot_style
set(w_ax(1:3),'Units',Plot_Units)
set(w_ax(1:3),'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2;
set(gcf,'DefaultLineLineWidth',Line_Width)
for n = 1:3
   set(w_ax(n),'Position',[Plot_X+(n-1)*Plot_Height/3,Plot_Y,Plot_Height/3,Plot_Height])
   axis(w_ax(n),[0 9 wmin wmax]);
   xlabel(w_ax(n),'P (bar)','Interpreter','LaTeX','FontSize',Label_Font_Size)
   title(w_ax(n),['Nozzle ' Nozzle{n}],'Interpreter','LaTeX','FontSize',Label_Font_Size)
   hl(n) = legend(hw(n,:),'Exp.','FDS','Location',w_leg_pos{n});
end
set(w_ax(2:3),'YTickLabel',[])
set(w_ax(1),'XTick',[0 2 4 6 8])
set(w_ax(2:3),'XTick',[2 4 6 8])
ylabel(w_ax(1),'Mean W-velocity (m/s)','Interpreter','LaTeX','FontSize',Label_Font_Size)

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Height Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Height Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_5_Validation_Guide/FIGURES/BRE_Spray/BRE_Spray_W');


% Attenuation plot

figure(hf(3))
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
Plot_X = 1.35*(Paper_Height-Plot_Height)/2;
Plot_Y = 1.25*(Paper_Height-Plot_Height)/2;
set(gca,'Position',[Plot_X,Plot_Y,Plot_Height,Plot_Height])
set(hf(3),'DefaultLineLineWidth',Line_Width)
xlabel('Exp. Attenuation (\%)','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('FDS Attenuation (\%)','Interpreter','LaTeX','FontSize',Label_Font_Size)
hl(n) = legend(ha,'Nozzle A','Nozzle B','Nozzle D','Location','SouthEast');

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Height Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Height Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_5_Validation_Guide/FIGURES/BRE_Spray/BRE_Spray_Attenuation');

close all
%EOF