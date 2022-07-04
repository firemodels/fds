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
pos = [19.8 321 426 508 554 604 628 682 723 834];
z_mod = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0 7.4]};
z_exp = {[0.3 1.1 2.0 2.6 3.2 3.7 4.1],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0]};
test  = {'612B','606A','615B'};

plot_style
Plot_Width      = 6.2;
Plot_Height     = 1.0;
Plot_X          = 0.2;
Plot_Y          = 0.2;
Paper_Width     = 6.5;
Paper_Height    = 1.3;

for k=1:3  % Experiment

   clear M E mod_data_indices exp_data_indices

   % Columns corresponding to velocity data at Loop positions: 214, 209, 208, 207, 307, 305, 304, 302, 301, 202
   mod_data_indices = {[9:15],[43:51],[61:69],[79:87],[97:105],[124:132],[151:159],[178:186],[196:204],[212:218]};
   exp_data_indices = {[73:79],[65:72],[57:64],[49:56],[41:48],[33:40],[25:32],[17:24],[9:16],[2:8]};

   M = importdata([outdir,'Test_',test{k},'_devc.csv'],',',2);
   E = importdata([expdir,'VP-',test{k},'.csv'],',',2);

   scale = -3;

   for i=1:19   % Time

      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*time(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*time(i),'nearest');

      for j=1:10   % Position

         clear u_mod u_exp
         u_mod = pos(j) + scale*M.data(mod_time_index,mod_data_indices{j});
         u_exp = pos(j) + scale*E.data(exp_time_index,exp_data_indices{j});
         if j==1 | j==10
            kk = 1;
         else
            kk = 2;
         end
         plot(u_mod,z_mod{kk},'k-') ; hold on
         plot(u_exp,z_exp{kk},'r-') ; hold on
         if j==1 | j==10
            plot([pos(j) pos(j)],[0 4.3],'k:') ; hold on
         else
            plot([pos(j) pos(j)],[0 8],'k:') ; hold on
         end
         plot([0 21.3 21.3],[4.3 4.3 8],'k-') ; hold on
         plot([854 832.7 832.7],[4.3 4.3 8],'k-') ; hold on

      end

      a = get(gca,'XTickLabel');  
      set(gca,'XTickLabel',a,'fontsize',6)
      set(gca,'TickLength',[0 0])
      xticks(pos)
      xticklabels({'214','209','208','207','307','305','304','302','301','202'})
      ylabel('Height (m)','FontSize',6,'Interpreter',Font_Interpreter)
      set(gca,'Units',Plot_Units)
      set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
      set(gca,'FontName',Font_Name)
      axis([0 854 0 8])
      text(30,7.2,['Test ',test{k},' Velocity Profiles'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      text(30,6.5,['Time: ',num2str(time(i)),' min'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)

      Git_Filename = [outdir,'Test_',test{k},'_git.txt'];
      addverstr(gca,Git_Filename,'linear',0.8,1.05,'Times','TeX',6)

      set(gcf,'Visible',Figure_Visibility);
      set(gcf,'Units',Paper_Units);
      set(gcf,'PaperUnits',Paper_Units);
      set(gcf,'PaperSize',[Paper_Width Paper_Height]);
      set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
      print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_U_',num2str(time(i))])

      hold off

   end


   % Temperature profiles

   clear E mod_data_indices exp_data_indices
   
   E = importdata([expdir,'TP-',test{k},'.csv'],',',2);
   
   % Columns corresponding to temperature data at Loop positions: 214, 209, 208, 207, 307, 305, 304, 302, 301, 202
   mod_data_indices = {[2:8],[34:42],[52:60],[70:78],[88:96],[115:123],[142:150],[169:177],[187:195],[205:211]};
   exp_data_indices = {[113:119],[89:96],[81:88],[73:80],[65:72],[49:56],[33:40],[17:24],[9:16],[2:8]};

   scale = 0.1;

   for i=1:19

      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*time(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*time(i),'nearest');
   
      for j=1:10
   
         clear u_mod u_exp
         u_mod = pos(j) + scale*M.data(mod_time_index,mod_data_indices{j});
         u_exp = pos(j) + scale*E.data(exp_time_index,exp_data_indices{j});
         if j==1 | j==10
            kk = 1;
         else
            kk = 2;
         end
         plot(u_mod,z_mod{kk},'k-') ; hold on
         plot(u_exp,z_exp{kk},'r-') ; hold on
         if j==1 | j==10
            plot([pos(j) pos(j)],[0 4.3],'k:') ; hold on
         else
            plot([pos(j) pos(j)],[0 8],'k:') ; hold on
         end
         plot([0 21.3 21.3],[4.3 4.3 8],'k-') ; hold on
         plot([854 832.7 832.7],[4.3 4.3 8],'k-') ; hold on

      end

      a = get(gca,'XTickLabel');  
      set(gca,'XTickLabel',a,'fontsize',6)
      set(gca,'TickLength',[0 0])
      xticks(pos)
      xticklabels({'214','209','208','207','307','305','304','302','301','202'})
      ylabel('Height (m)','FontSize',6,'Interpreter',Font_Interpreter)
      set(gca,'Units',Plot_Units)
      set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
      set(gca,'FontName',Font_Name)
      axis([0 854 0 8])
      text(30,7.2,['Test ',test{k},' Temperature Profiles'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)
      text(30,6.5,['Time: ',num2str(time(i)),' min'],'Fontname',Font_Name,'FontSize',6,'Interpreter',Font_Interpreter)

      Git_Filename = [outdir,'Test_',test{k},'_git.txt'];
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

M = importdata([outdir,'Cold_Flow_Series_1_devc.csv'],',',2);
E = importdata([expdir,'Cold_Flow_Series_1.csv'],',',2);

for j=1:15
   mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),300*j,'nearest');
end

figure

plot(E.data(:,1),E.data(:,3),'ko-') ; hold on
plot(M.data(mod_time_index,1)/300,M.data(mod_time_index,2),'ko--') ; hold on

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 16 0 500])
legend_handle=legend('Measured','FDS','Location','SouthEast');
set(legend_handle,'Fontsize',Key_Font_Size);
xlabel('Number of Fans','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Volume Flow (m³/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
Git_Filename = [outdir,'Cold_Flow_Series_1_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Cold_Flow_Series_1_Volume_Flow'])

display('Memorial_Tunnel completed successfully')
