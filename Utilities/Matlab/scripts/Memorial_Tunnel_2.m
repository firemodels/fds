% McGrattan
% 11-17-2023
% Memorial_Tunnel_2.m
%
% This script creates several different kinds of contour and scatter plots for the Memorial Tunnel simulations with the ceiling in place.

clear all
close all

outdir = '../../../out/Memorial_Tunnel/';
expdir = '../../../exp/Memorial_Tunnel/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/';

times = [1];

% Loop  214  213   211   209   208   207   307   306   305   205   304   303   302   301   202
pos  = [19.8 125.6 230.7 340.8 446.2 528.2 573.3 586.1 604.1 615.4 627.6 645.0 681.5 723.3 833.9];
vpos = [19.8             340.8 446.2 528.2 573.3       604.1       627.6       681.5 723.3 833.9];
hgt_mod = {[0.3 1.1 1.8 2.4 3.0 3.5 4.0],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0 7.4]};
hgt_exp = {[0.3 1.1 1.8 2.4 3.0 3.5 4.0],[0.3 1.2 2.4 3.7 4.8 5.7 6.5 7.0]};

test  = { '101CR', '102', '102R1', '102R', '103', '104', '105', '106', '107', '108',...
          '109', '110', '111', '112A', '113A', '115A', '126B', '126BR1', '128B', '202',...
          '203', '205', '207A', '208A', '210', '212', '214A', '215A', '216A', '217A',...
          '218B', '223', '226', '227A', '229', '230', '231', '233', '235', '236',...
          '238A', '239', '244B', '245B', '246B', '247B', '248B', '249B', '250B', '251B',...
          '252B', '301A', '302A', '303A', '305A', '306A', '309A', '312A', '313A', '314',...
          '315A', '316', '317A', '318A', '319A', '320A', '321A', '338B', '339B', '340B',...
          '341B', '342B', '343B', '344B', '345B', '346B', '401A', '403A', '404A', '407B',...
          '408B'};

sequence = [ 1 1 1 1 1 1 1 1 1 1 ...
             1 1 1 1 1 1 1 1 1 3 ...
             3 3 3 3 3 3 3 3 3 3 ...
             3 4 5 5 6 6 6 6 6 6 ...
             6 6 5 6 6 3 3 3 6 6 ...
             3 8 8 8 8 8 8 8 8 9 ...
             101 9 102 101 102 102 9 8 8 8 ...
             8 8 8 8 8 8 13 13 14 13 ...
             13 ];

levels = [50 100 200 400 600 800];
single_level = [50 50];
% T Loops:           214     , 213     , 211    , 209   , 208   , 207    , 307     , 306     , 305     , 205     , 304     , 303     , 302     , 301     , 202
mod_data_indices = {[2:8]    ,[16:22]  ,[23:29] ,[30:36],[44:50],[58:64] ,[72:78]  ,[86:92]  ,[93:99]  ,[107:113],[114:120],[128:134],[135:141],[149:155],[163:169]};
exp_data_indices = {[100:106],[93:99]  ,[86:92] ,[79:85],[72:78],[65:71] ,[58:64]  ,[51:57]  ,[44:50]  ,[37:43]  ,[30:36]  ,[23:29]  ,[16:22]  ,[9:15]   ,[2:8]};

v_loop =           {'214'  ,'209'    ,'208'   ,'207'   ,'307'    ,'305'     ,'304'     ,'302'     ,'301'     ,'202'};
mod_vel_indices = {[9:15]  ,[37:43]  ,[51:57] ,[65:71] ,[79:85]  ,[100:106] ,[121:127] ,[142:148] ,[156:162] ,[170:176]};
exp_vel_indices = {[65:71] ,[58:64]  ,[51:57] ,[44:50] ,[37:43]  ,[30:36]   ,[23:29]   ,[16:22]   ,[9:15]    ,[2:8]};

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

for k=1:81 % Experiments

   clear M E EV H
   
   try
      M  = importdata([outdir,'Test_',test{k},'_cat_devc.csv'],',',2);
      E  = importdata([expdir,'TP',test{k},'.csv'],',',2);
      EV = importdata([expdir,'QP',test{k},'.csv'],',',2);
      H  = importdata([expdir,'HR',test{k},'.csv'],',',2);
   catch
      continue
   end

   if M.data(end,1)<200 ; continue ; end
   
   switch sequence(k)
      case 1      ; ventilation_type = 'Full Transverse Ventilation';
      case 3      ; ventilation_type = 'Partial Transverse Exhaust Ventilation';
      case {4,5}  ; ventilation_type = 'Partial Transverse Exhaust Ventilation';
      case 6      ; ventilation_type = 'Two-Zone Partial Transverse Ventilation';
      case 8      ; ventilation_type = 'Partial Transverse Ventilation with Single Point Extraction';
      case 9      ; ventilation_type = 'Point Supply Operation';
      case 10     ; ventilation_type = 'Point Exhaust Operation';
      case {13,14}; ventilation_type = 'Partial Transverse Ventilation with Oversize Exhaust Ports';
   end

   for i=1:length(times) % Times
       
      clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp
       
      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*times(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*times(i),'nearest');
      hrr_time_index = interp1(H.data(:,1),1:length(H.data(:,1)),60*times(i),'nearest');
      exp_VF_time_index = interp1(EV.data(:,1),1:length(EV.data(:,1)),60*times(i),'nearest');
   
      [X_mod,Y_mod] = meshgrid(pos(1:15),hgt_mod{1});
      [X_exp,Y_exp] = meshgrid(pos(1:15),hgt_exp{1});
      for kk=1:15 % Loops
         Z_mod(:,kk) = M.data(mod_time_index,mod_data_indices{kk});
         Z_exp(:,kk) = E.data(exp_time_index,flip(exp_data_indices{kk}));
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
      axis([0 854 0 4.4])
      text(10,4.0,['Test ',test{k}],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,3.5,['Time: ',num2str(times(i)),' min'],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,3.0,['HRR: ',num2str(H.data(hrr_time_index,2)/1000.,'%.1f'),' MW'],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      text(10,2.5,'FDS red; Exp black','Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
      
      Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
      addverstr(gca,Git_Filename,'linear',0.75,1.07,'Times','TeX',7)

      set(gcf,'Units',Paper_Units);
      set(gcf,'PaperUnits',Paper_Units);
      set(gcf,'PaperSize',[Paper_Width Paper_Height+0.2]);
      set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
      print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_T_',num2str(times(i))])

      hold off

   end

   % For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

   clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp

   [X_mod,Y_mod] = meshgrid(pos(1:15),M.data(:,1)/60);
   [X_exp,Y_exp] = meshgrid(pos(1:15),E.data(:,1)/60);
   for kk=1:length(M.data(:,1))
      for ii=1:15
         Z_mod(kk,ii) = M.data(kk,mod_data_indices{ii}(7));
      end
   end

   for kk=1:length(E.data(:,1))
      for ii=1:15
         Z_exp(kk,ii) = E.data(kk,exp_data_indices{ii}(1));
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
   xticks([0 100 200 300 400 500 600 700 800 854]);
   xticklabels({'0','100','200','300','400','500','600','700','800','854'})
   xlabel('Tunnel Length (m)','FontSize',7,'Interpreter',Font_Interpreter)
   ylabel('Time (min)','FontSize',7,'Interpreter',Font_Interpreter)
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X_2 (Plot_Y_2+0.1) Plot_Width_2 Plot_Height_2])
   set(gca,'FontName',Font_Name)
   axis([0 854 0 30])
   text(10,26,['Test ',test{k}],'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
   text(10,22,ventilation_type,'Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)
   text(10,18,'FDS red; Exp black','Fontname',Font_Name,'FontSize',7,'Interpreter',Font_Interpreter)

   Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
   addverstr(gca,Git_Filename,'linear',0.75,1.07,'Times','TeX',7)

   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height+0.2]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_tvT'])

   hold off

   % Make velocity profile plots

   clear M E EV H
   
   M  = importdata([outdir,'Test_',test{k},'_cat_devc.csv'],',',2);
   E  = importdata([expdir,'VP',test{k},'.csv'],',',2);
   
   for i=1:length(times) % Times
       
      clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp
       
      mod_time_index = interp1(M.data(:,1),1:length(M.data(:,1)),60*times(i),'nearest');
      exp_time_index = interp1(E.data(:,1),1:length(E.data(:,1)),60*times(i),'nearest');
   
      for kk=1:10 % Loops
         v_mod = M.data(mod_time_index,mod_vel_indices{kk});
         v_exp = E.data(exp_time_index,flip(exp_vel_indices{kk}));
         plot(vpos(kk)-5*v_mod,hgt_mod{1},'r-') ; hold on
         plot(vpos(kk)-5*v_exp,hgt_exp{1},'k-') ; hold on
         plot([vpos(kk) vpos(kk)],[0 4.4],'k-') ; hold on
      end

      a = get(gca,'XTickLabel');  
      set(gca,'XTickLabel',a,'fontsize',7)
      set(gca,'TickLength',[0 0])
      xticks(vpos);
      xticklabels(v_loop)
      xlabel('Velocity Measurement Locations','FontSize',7,'Interpreter',Font_Interpreter)
      ylabel('Height (m)','FontSize',7,'Interpreter',Font_Interpreter)
      set(gca,'Units',Plot_Units)
      set(gca,'Position',[Plot_X_2 (Plot_Y_2+0.1) Plot_Width_2 Plot_Height_2])
      set(gca,'FontName',Font_Name)
      axis([0 854 0 4.4])
      
      Git_Filename = [outdir,'Test_',test{k},'_cat_git.txt'];
      addverstr(gca,Git_Filename,'linear',0.75,1.07,'Times','TeX',7)

      set(gcf,'Units',Paper_Units);
      set(gcf,'PaperUnits',Paper_Units);
      set(gcf,'PaperSize',[Paper_Width Paper_Height+0.2]);
      set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
      print(gcf,'-dpdf',[pltdir,'Test_',test{k},'_V_',num2str(times(i))])

      hold off

   end

end  % Experiment loop

display('Memorial_Tunnel_2 completed successfully')
