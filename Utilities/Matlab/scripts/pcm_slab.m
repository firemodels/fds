% Hostikka
% 11-15-2010
% pcm_clab.m
%
% This script reads and post-processes the output from
% a phase-change verification test pcm_slab.fds
% and plots comparisons against the analytical results
%


close all
clear all

% Calculate analytical result

t = 0:60:3.6E4;
x = [0.01 0.05 0.1];
lambda = 0.189;
alpha = 2/(1000*2000);
T_f = 0;
T_0 = -15;
T_i = 15;
%
xf = lambda*2*(alpha*t).^(1/2);
%
nt = length(t);
for n = 1:nt
   for i = 1:length(x)
      l_n = x(i)/(2*sqrt(alpha*t(n)));
      if (x(i) < xf(n))
         ER = erf(l_n)/erf(lambda);
         T(n,i) = T_0 + (T_f-T_0)*ER;
      else
         ER = erfc(l_n)/erfc(lambda);
         T(n,i) = T_i - (T_i-T_f)*ER;
      end
   end
end

addpath('../../Verification/Pyrolysis')

skip_case = 0;
if ~exist('pcm_slab_prof_01.csv')
    display('Error: File pcm_slab_prof_01.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('pcm_slab_devc.csv')
    display('Error: File pcm_slab_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if skip_case
    return
end

Iprof = csvread('pcm_slab_prof_01.csv',2);
M_fds = csvread('pcm_slab_devc.csv',2);
%
Ice_density = 990;
[nt_fds,nc] = size(Iprof);
Nw = Iprof(1,2);
%
x_fds = Iprof(1,3:Nw+2);
t_fds = Iprof(:,1);
Iprof=Iprof(:,[1 Nw+3:nc]);
%
[nt_fds,nc] = size(Iprof);
for n = 1:nt_fds
   i1 = find(Iprof(n,2:nc)>0,1);
   i2 = find(Iprof(n,2:nc)>0,1,'last');
   if (isempty(i1))
      dmax(n) = 0;
   else
      dmax(n) = x_fds(i2)-x_fds(i1);
   end
   i1 = find(Iprof(n,2:nc)>Ice_density,1);
   i2 = find(Iprof(n,2:nc)>Ice_density,1,'last');
   if (isempty(i1))
      dmin(n) = 0;
   else
      dmin(n) = x_fds(i2)-x_fds(i1);
   end
end

for i_plot=1:2
    
    close all
    
    plot_style
    set(gcf,'DefaultLineLineWidth',Line_Width)
    set(gca,'FontName',Font_Name)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

    if i_plot==1
       h=plot(t,xf,'k-',t_fds,dmin,'k--',t_fds,dmax,'r--');
       xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
       ylabel('Phase interface, {\itx}_f (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
       legend('Analytical','FDS min','FDS max','Location','SouthEast')
    else
       h=plot(t,T(:,1),'k-',t,T(:,2),'r-',t,T(:,3),'g-'); hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,2),'kd');hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,3),'ro');hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,4),'gs');
       xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
       ylabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)       
       legend('Analytical 1 cm','Analytical 5 cm','Analytical 10 cm','FDS 1 cm','FDS 5 cm','FDS 10 cm')
    end
    
    % Plot attributes
    
    set(gca,'FontName',Font_Name)
    
    % add SVN if file is available
    
    SVN_Filename = ['pcm_slab_svn.txt'];
    if exist(SVN_Filename,'file')
        SVN = importdata(SVN_Filename);
        x_lim = get(gca,'XLim');
        y_lim = get(gca,'YLim');
        X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
        Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
        text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
            'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
    end
    
    % Create the PDF files
    
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
    if i_plot==1
        print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/pcm_slab_xf')
    else
        print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/pcm_slab_T')
    end
    
end
