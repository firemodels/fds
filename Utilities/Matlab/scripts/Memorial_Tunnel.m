% McGrattan
% 6-23-2022
% Memorial_Tunnel.m
%
% This script creates vertical profiles of velocity and temperature for the Memorial Tunnel simulations.

clear all
close all

outdir = '../../../out/Memorial_Tunnel/';
expdir = '../../../exp/Memorial_Tunnel/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/';

time = [0 1 2 3 4 5 6 8 10 12 14 16 18 20 22 24 26 28 30];
pos = [19.8 125.6 230.7 340.8 446.2 528.2 573.3 586.1 604.1 615.4 627.6 645.0 681.5 723.3 833.9];
hgt_mod = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0 7.4]};
hgt_exp = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0]};
test  = {'501','502','605','606A','607','608','610','611','612B','615B','617A','618A','621A','622B','623B','624B','625B'};
%levels = [21 38 60 93 149 204 260 316 371 427 482 538 593 649 704 760 816];  % Original values 
levels = [50 100 200 400 600 800];
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

for k=1:17 % Experiments

   clear M E
   
   M = importdata([outdir,'Test_',test{k},'_cat_devc.csv'],',',2);
   E = importdata([expdir,'TP-',test{k},'.csv'],',',2);
   H = importdata([expdir,'HRR-',test{k},'.csv'],',',2);
   
   for i=1:18 % Times
       
      clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp
       
      if k==13 & i>17 ; break ; end
      if k==17 & i>16 ; break ; end

      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*time(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*time(i),'nearest');
      hrr_time_index = interp1(H.data(:,1),1:length(H.data(:,1)),60*time(i),'nearest');
   
      [X_mod,Y_mod] = meshgrid(pos(2:14),hgt_mod{2});
      [X_exp,Y_exp] = meshgrid(pos(2:14),hgt_exp{2});
      for kk=2:14 % Loops
         Z_mod(:,kk-1) = M.data(mod_time_index,mod_data_indices{kk});
         Z_exp(:,kk-1) = E.data(exp_time_index,exp_data_indices{kk});
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

      [C_mod,h_mod] = contour(X_mod_interp,Y_mod_interp,Z_mod_interp,levels,'k-') ; hold on
      [C_exp,h_exp] = contour(X_exp_interp,Y_exp_interp,Z_exp_interp,levels,'r-') ; hold on
      clabel(C_mod,h_mod,'FontSize',3,'Color','black','LabelSpacing',300)
      clabel(C_exp,h_exp,'FontSize',3,'Color','red','LabelSpacing',300)

      a = get(gca,'XTickLabel');  
      set(gca,'XTickLabel',a,'fontsize',4)
      set(gca,'TickLength',[0 0])
      xticks(pos)
      xticklabels({'214','213','211','209','208','207','307','306','305','F','304','303','302','301','202'})
      ylabel('Height (m)','FontSize',5,'Interpreter',Font_Interpreter)
      set(gca,'Units',Plot_Units)
      set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
      set(gca,'FontName',Font_Name)
      axis([0 854 0 8])
      text(10,7.2,['Test ',test{k}],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      text(10,6.3,['Time: ',num2str(time(i)),' min'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      text(10,5.4,['HRR: ',num2str(H.data(hrr_time_index,2)/1000.,'%.1f'),' MW'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      text(10,4.5,'FDS in black; Exp in red','Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      
      Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
      addverstr(gca,Git_Filename,'linear',0.8,1.05,'Times','TeX',6)

      %set(gcf,'Visible',Figure_Visibility);
      set(gcf,'Units',Paper_Units);
      set(gcf,'PaperUnits',Paper_Units);
      set(gcf,'PaperSize',[Paper_Width Paper_Height]);
      set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
      print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_T_',num2str(time(i))])

      hold off

   end

end

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

plot([3 3 9 15 15],[169.4 164.7 292.6 372.4 379.9],'r^') ; hold on
%plot(E.data(:,1),E.data(:,3),'ko-') ; hold on
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
