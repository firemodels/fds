% McGrattan
% 7-29-2022
% Memorial_Tunnel.m
%
% This script creates several different kinds of contour and scatter plots for the Memorial Tunnel simulations.

clear all
close all

outdir = '../../../out/Memorial_Tunnel/';
expdir = '../../../exp/Memorial_Tunnel/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/';

times = {[ 5 10 16 20 26 30],...  % 501
         [ 5 10 16 20 26 30],...  % 502
         [12 16 22 26 28 30],...  % 605
         [ 4  8 16 20 24 26],...  % 606A
         [ 8 10 14 22 26 28],...  % 607
         [ 1  2  4  6 12 20],...  % 608
         [ 6  8 12 22 24 26],...  % 610
         [ 2  6  8 16 20 26],...  % 611
         [ 2  6 10 14 18 22],...  % 612B
         [ 2  4  6 10 18 28],...  % 615B
         [16 18 20 22 28 30],...  % 617A
         [ 2  5  6 12 24 26],...  % 618A
         [ 1  2  3 12 14 24],...  % 621A
         [ 3  4  5  8 10 14],...  % 622B
         [ 2  4  6 10 20 28],...  % 623B
         [ 1  2  3  4  6  8],...  % 624B
         [ 6  8 10 18 20 21]};    % 625B
pos = [19.8 125.6 230.7 340.8 446.2 528.2 573.3 586.1 604.1 615.4 627.6 645.0 681.5 723.3 833.9];
hgt_mod = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0 7.4]};
hgt_exp = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0]};
test  = {'501','502','605','606A','607','608','610','611','612B','615B','617A','618A','621A','622B','623B','624B','625B'};
levels = [50 100 200 400 600 800];
single_level = [50 50];
% Loops:             214     , 213     , 211    , 209   , 208   , 207    , 307     , 306     , 305     , 205     , 304     , 303     , 302     , 301     , 202
mod_data_indices = {[2:8]    ,[16:24]  ,[25:33] ,[34:42],[52:60],[70:78] ,[88:96]  ,[106:114],[115:123],[133:141],[142:150],[160:168],[169:177],[187:195],[205:211]};
exp_data_indices = {[113:119],[105:112],[97:104],[89:96],[81:88],[73:80] ,[65:72]  ,[57:64]  ,[49:56]  ,[41:48]  ,[33:40]  ,[25:32]  ,[17:24]  ,[9:16]   ,[2:8]};

plot_style
Plot_Width      = 6.2;
Plot_Height     = 1.0;
Plot_X          = 0.2;
Plot_Y          = 0.2;
Paper_Width     = 6.5;
Paper_Height    = 1.3;

Plot_Width_2    = 6.15;
Plot_Height_2   = 1.0;
Plot_X_2        = 0.3;
Plot_Y_2        = 0.1;

% Make contour plots of centerline tunnel temperature for 17 experiments and 19 individual times.
% Also, make a contour plot for each experiment showing a single temperature contour as a function of position and time.

for k=1:17 % Experiments

   clear M E EV H
   
   M  = importdata([outdir,'Test_',test{k},'_cat_devc.csv'],',',2);
   E  = importdata([expdir,'TP',test{k},'.csv'],',',2);
   EV = importdata([expdir,'QP',test{k},'.csv'],',',2);
   H  = importdata([expdir,'HRR',test{k},'.csv'],',',2);
   
   for i=1:length(times{k}) % Times
       
      clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp
       
      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*times{k}(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*times{k}(i),'nearest');
      hrr_time_index = interp1(H.data(:,1),1:length(H.data(:,1)),60*times{k}(i),'nearest');
      exp_VF_time_index = interp1(EV.data(:,1),1:length(EV.data(:,1)),60*times{k}(i),'nearest');
   
      [X_mod,Y_mod] = meshgrid(pos(2:14),hgt_mod{2});
      [X_exp,Y_exp] = meshgrid(pos(2:14),hgt_exp{2});
      for kk=2:14 % Loops
         Z_mod(:,kk-1) = M.data(mod_time_index,mod_data_indices{kk});
         Z_exp(:,kk-1) = E.data(exp_time_index,flip(exp_data_indices{kk}));
      end

      newpoints = 100;
      [X_mod_interp,Y_mod_interp] = meshgrid(...
            linspace(min(min(X_mod,[],2)),max(max(X_mod,[],2)),newpoints ),...
            linspace(min(min(Y_mod,[],1)),max(max(Y_mod,[],1)),newpoints )...
          );
      Z_mod_interp = interp2(X_mod,Y_mod,Z_mod,X_mod_interp,Y_mod_interp,'makima');
      
      [X_exp_interp,Y_exp_interp] = meshgrid(...
            linspace(min(min(X_exp,[],2)),max(max(X_exp,[],2)),newpoints ),...
            linspace(min(min(Y_exp,[],1)),max(max(Y_exp,[],1)),newpoints )...
          );
      Z_exp_interp = interp2(X_exp,Y_exp,Z_exp,X_exp_interp,Y_exp_interp,'makima');

      reset(gca)
      reset(gcf)

      [C_mod,h_mod] = contour(X_mod_interp,Y_mod_interp,Z_mod_interp,levels,'r-') ; hold on
      [C_exp,h_exp] = contour(X_exp_interp,Y_exp_interp,Z_exp_interp,levels,'k-') ; hold on
      clabel(C_mod,h_mod,'FontSize',3,'Color','red'  ,'LabelSpacing',300)
      clabel(C_exp,h_exp,'FontSize',3,'Color','black','LabelSpacing',300)

      a = get(gca,'XTickLabel');  
      set(gca,'XTickLabel',a,'fontsize',7)
      set(gca,'TickLength',[0 0])
      xticks([0 100 200 300 400 500 600 700 800 854]);
      xticklabels({'0','100','200','300','400','500','600','700','800','854'})
      xlabel('Tunnel Length (m)','FontSize',7,'Interpreter',Font_Interpreter)
      ylabel('Height (m)','FontSize',7,'Interpreter',Font_Interpreter)
      set(gca,'Units',Plot_Units)
      set(gca,'Position',[Plot_X_2 (Plot_Y_2+0.1) Plot_Width_2 Plot_Height_2])
      set(gca,'FontName',Font_Name)
      axis([0 854 0 8])
      text(10,7.2,['Test ',test{k}],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,6.3,['Time: ',num2str(times{k}(i)),' min'],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,5.4,['HRR: ',num2str(H.data(hrr_time_index,2)/1000.,'%.1f'),' MW'],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,4.5,'FDS red; Exp black','Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      
      Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
      addverstr(gca,Git_Filename,'linear',0.75,1.07,'Times','TeX',7)

      set(gcf,'Units',Paper_Units);
      set(gcf,'PaperUnits',Paper_Units);
      set(gcf,'PaperSize',[Paper_Width Paper_Height+0.2]);
      set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
      print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_T_',num2str(times{k}(i))])

      hold off

   end

   % For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

   clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp

   [X_mod,Y_mod] = meshgrid(pos(2:14),M.data(:,1)/60);
   [X_exp,Y_exp] = meshgrid(pos(2:14),E.data(:,1)/60);
   for kk=1:length(M.data(:,1))
      for ii=2:14
         Z_mod(kk,ii-1) = M.data(kk,mod_data_indices{ii}(8));
      end
   end

   for kk=1:length(E.data(:,1))
      for ii=2:14
         Z_exp(kk,ii-1) = E.data(kk,exp_data_indices{ii}(1));
      end
   end

   newpoints = 100;
   [X_mod_interp,Y_mod_interp] = meshgrid(...
         linspace(min(min(X_mod,[],2)),max(max(X_mod,[],2)),newpoints ),...
         linspace(min(min(Y_mod,[],1)),max(max(Y_mod,[],1)),newpoints )...
       );
   Z_mod_interp = interp2(X_mod,Y_mod,Z_mod,X_mod_interp,Y_mod_interp,'makima');

   [X_exp_interp,Y_exp_interp] = meshgrid(...
         linspace(min(min(X_exp,[],2)),max(max(X_exp,[],2)),newpoints ),...
         linspace(min(min(Y_exp,[],1)),max(max(Y_exp,[],1)),newpoints )...
       );
   Z_exp_interp = interp2(X_exp,Y_exp,Z_exp,X_exp_interp,Y_exp_interp,'makima');

   reset(gca)
   reset(gcf)

   [C_mod,h_mod] = contour(X_mod_interp,Y_mod_interp,Z_mod_interp,single_level,'r-') ; hold on
   [C_exp,h_exp] = contour(X_exp_interp,Y_exp_interp,Z_exp_interp,single_level,'k-') ; hold on
   clabel(C_mod,h_mod,'FontSize',3,'Color','red'  ,'LabelSpacing',300)
   clabel(C_exp,h_exp,'FontSize',3,'Color','black','LabelSpacing',300)
   plot([615.4 615.4],[0 30],'k--')

   a = get(gca,'XTickLabel');
   set(gca,'XTickLabel',a,'fontsize',7)
   set(gca,'TickLength',[0 0])
   xticks(pos)
  %xticklabels({'214','213','211','209','208','207','307','306','305','F','304','303','302','301','202'})
   xticklabels({'','','','','','','','','','','','','','',''})
   ylabel('Time (min)','FontSize',7,'Interpreter',Font_Interpreter)
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X_2 Plot_Y_2 Plot_Width_2 Plot_Height_2])
   set(gca,'FontName',Font_Name)
   axis([0 854 0 30])
   text(10,26,['Test ',test{k}],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
   text(10,22,'FDS red; Exp black','Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)

   Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
   addverstr(gca,Git_Filename,'linear',0.75,1.07,'Times','TeX',7)

   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_tvT'])

   hold off

end  % Experiment loop


% Collect critical velocity information

cv_time = {[10 12 14 16 18 20],...  % 501
           [10 12 14 16 18 20],...  % 502
           [25 28],...            % 605
           [24],...  % 606A
           [2 10],...  % 607
           [ ],...  % 608
           [2 8 20 25],...  % 610 
           [10 16],...  % 611 
           [15 20],...  % 612B
           [13 20],...  % 615B
           [14 17 27 30],...  % 617A
           [20 26],...  % 618A
           [8 12],...  % 621A
           [3 6],...  % 622B
           [ ],...  % 623B
           [10 14],...  % 624B
           [ ]};    % 625B

js = 0;
jn = 0;
ks = 0;
kn = 0;

for k=3:17 % Experiments

   clear M E EV H

   M  = importdata([outdir,'Test_',test{k},'_cat_devc.csv'],',',2);
   E  = importdata([expdir,'TP',test{k},'.csv'],',',2);
   EV = importdata([expdir,'QP',test{k},'.csv'],',',2);
   H  = importdata([expdir,'HRR',test{k},'.csv'],',',2);
   HM = importdata([outdir,'Test_',test{k},'_cat_hrr.csv'],',',2);

   for i=1:length(cv_time{k}) % Times

      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*cv_time{k}(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*cv_time{k}(i),'nearest');
      hrr_time_index = interp1(H.data(:,1),1:length(H.data(:,1)),60*cv_time{k}(i),'nearest');
      hrr_mod_time_index = interp1(HM.data(:,1),1:length(HM.data(:,1)),60*cv_time{k}(i),'nearest');
      exp_VF_time_index = interp1(EV.data(:,1),1:length(EV.data(:,1)),60*cv_time{k}(i),'nearest');

      if E.data(exp_time_index,64)<30
         js = js+1;
         V_crit_exp_sup(js) = -EV.data(exp_VF_time_index,2)/60.4;
         HRR_crit_exp_sup(js) = H.data(hrr_time_index,2)/1000;
      else
         jn = jn+1;
         V_crit_exp_nosup(jn) = -EV.data(exp_VF_time_index,2)/60.4;
         HRR_crit_exp_nosup(jn) = H.data(hrr_time_index,2)/1000;
      end

      if M.data(mod_time_index,113)<30
         ks = ks+1;
         V_crit_mod_sup(ks) = -M.data(mod_time_index,219)/60.4;
         HRR_crit_mod_sup(ks) = HM.data(hrr_mod_time_index,2)/1000;
      else
         kn = kn+1;
         V_crit_mod_nosup(kn) = -M.data(mod_time_index,219)/60.4;
         HRR_crit_mod_nosup(kn) = HM.data(hrr_mod_time_index,2)/1000;
      end

   end

end

figure

reset(gca)
reset(gcf)
plot_style

semilogx(HRR_crit_exp_sup,V_crit_exp_sup,'kd') ; hold on
semilogx(HRR_crit_mod_sup,V_crit_mod_sup,'rd') ; hold on
h=semilogx(HRR_crit_exp_nosup,V_crit_exp_nosup,'ks') ; hold on
set(h,'MarkerFaceColor',get(h,'Color'));
h2=semilogx(HRR_crit_mod_nosup,V_crit_mod_nosup,'rs') ; hold on
set(h2,'MarkerFaceColor',get(h2,'Color'));

semilogx([8.5 105],[1.9 3.5],'k-') ; hold on

xtick = [10 50 100];
xticks(xtick);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([5 200 0 4])
legend_handle=legend('Backlayering Controlled (Exp)','Backlayering Controlled (FDS)',...
                     'Backlayering Not Controlled (Exp)','Backlayering Not Controlled (FDS)','Theoretical Critical Velocity','Location','SouthEast');
set(legend_handle,'Fontsize',Key_Font_Size);
xlabel('Heat Release Rate (MW)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Velocity (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
Git_Filename = [outdir,'Cold_Flow_Series_1_cat_git.txt'];
addverstr(gca,Git_Filename,'semilogx')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Critical_Velocity'])


% Process Cold Flow case

clear all

outdir = '../../../out/Memorial_Tunnel/';
expdir = '../../../exp/Memorial_Tunnel/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/';

plot_style

M = importdata([outdir,'Cold_Flow_Series_1_cat_devc.csv'],',',2);
E = importdata([expdir,'Cold_Flow_Series_1.csv'],',',2);

for j=1:15
   mod_time_index(j) = interp1(M.data(:,1),1:length(M.data(:,1)),300*j,'nearest');
end

figure

plot([3 3 9 15 15],[169.4 164.7 292.6 372.4 379.9],'k^') ; hold on
plot(M.data(mod_time_index,1)/300,M.data(mod_time_index,2),'ko-') ; hold on

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 16 0 400])
legend_handle=legend('Measured','FDS','Location','SouthEast');
set(legend_handle,'Fontsize',Key_Font_Size);
xlabel('Number of Fans','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Volume Flow (m³/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
xticks([1 3 5 7 9 11 13 15]);

% add Git revision if file is available
Git_Filename = [outdir,'Cold_Flow_Series_1_cat_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Cold_Flow_Series_1_Volume_Flow'])

display('Memorial_Tunnel completed successfully')
