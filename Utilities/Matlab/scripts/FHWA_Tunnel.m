% McGrattan
% 9-27-2022
% FHWA_Tunnel.m
%
% This script creates several different kinds of contour and scatter plots for the FHWA Tunnel simulations.

clear all
close all

outdir = '../../../out/FHWA_Tunnel/';
expdir = '../../../exp/FHWA_Tunnel/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FHWA_Tunnel/';

pos = [1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.5 9.5 10.5 11.5];
test  = {'IFAB-07','IFAB-08','IFAB-09','IFAB-10','IFAB-11','IFAB-13','IFAB-14','IFAB-15','IFAB-19','IFAB-22','IFAB-24'};
test2 = {'Test 7','Test 8','Test 9','Test 10','Test 11','Test 13','Test 14','Test 15','Test 19','Test 22','Test 24'};
single_level = [50 50];
setpoint = [10000 400 399 338 240 322 390 420 360 10000 10000];

plot_style

%Plot_Width    = 6.2;
%Plot_Height   = 1.0;
%Plot_X        = 0.2;
%Plot_Y        = 0.2;
%Paper_Width   = 6.5;
%Paper_Height  = 1.3;

for k=1:11 % Experiments

   n_res = 1;
%  if k==8 
%     n_res = 2;
%  end

   for jj=1:n_res

      clear M E

      M  = importdata([outdir,test{k},'_cat_devc.csv'],',',2);
      if jj==1
         M  = importdata([outdir,test{k},'_cat_devc.csv'],',',2);
         E  = importdata([expdir,test{k},'_avg.csv'],',',2);
      elseif jj==2
         M  = importdata([outdir,test{k},'_fine_cat_devc.csv'],',',2);
         E  = importdata([expdir,test{k},'_avg.csv'],',',2);
      end
   
      % For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

      clear X_mod Y_mod Z_mod X_exp Y_exp Z_exp

      [X_mod,Y_mod] = meshgrid(pos(1:16),M.data(:,1)/60);
      [X_exp,Y_exp] = meshgrid(pos(1:16),E.data(:,1)/60);
      for kk=1:length(M.data(:,1))
         for ii=1:16
            Z_mod(kk,ii) = M.data(kk,ii+1);
         end
      end
   
      for kk=1:length(E.data(:,1))
         for ii=1:16
            Z_exp(kk,ii) = E.data(kk,ii+4);
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
   
      if jj==1
         reset(gca)
         reset(gcf)
         mod_symbol = 'r-';
      else
         mod_symbol = 'r--';
      end
   
      [C_mod,h_mod] = contour(X_mod_interp,Y_mod_interp,Z_mod_interp,single_level,mod_symbol) ; hold on
      clabel(C_mod,h_mod,'FontSize',3,'Color','red'  ,'LabelSpacing',300)
      [C_exp,h_exp] = contour(X_exp_interp,Y_exp_interp,Z_exp_interp,single_level,'k-') ; hold on
      clabel(C_exp,h_exp,'FontSize',3,'Color','black','LabelSpacing',300)
   
   end  % grid resolution cases

   plot([5.5 5.5],[0 15],'k--')
   plot([0.0 15.],[setpoint(k)/60 setpoint(k)/60],'k:')

  %xticks(pos)
  %xticklabels({'1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.5','9.5','10.5','11.5'})
   ax = gca;
   ax.XAxis.FontSize = 16;
   ax.YAxis.FontSize = 16;
   xlabel('Position (m)','FontSize',16,'Interpreter',Font_Interpreter)
   ylabel('Time (min)','FontSize',16,'Interpreter',Font_Interpreter)
   set(gca,'Units',Plot_Units)
   set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
   set(gca,'FontName',Font_Name)
   axis([0 15 0 15])
   text(0.5,4,[test2{k}],'Fontname',Font_Name,'FontSize',10,'Interpreter',Font_Interpreter)
   text(0.5,3,'FDS red; Exp black','Fontname',Font_Name,'FontSize',10,'Interpreter',Font_Interpreter)

   Git_Filename = [outdir,test{k},'_cat_git.txt'];
   addverstr(gca,Git_Filename,'linear',0.6,1.05,'Times','TeX',10)

   set(gcf,'Units',Paper_Units);
   set(gcf,'PaperUnits',Paper_Units);
   set(gcf,'PaperSize',[Paper_Width Paper_Height]);
   set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
   print(gcf,'-dpdf',[pltdir,test{k},'_tvT'])

   hold off

end  % Experiment loop


display('FHWA_Tunnel completed successfully')
