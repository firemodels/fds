% Hostikka
% 11-15-2010
% pcm_clab.m
%
% This script reads and post-processes the output from
% a phase-change verification test pcm_slab.fds
% and plots comparisons against the analytical results

close all
clear all

plot_style

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

    if i_plot==1
       h=plot(t,xf,'k-',t_fds,dmin,'b--',t_fds,dmax,'r--');
       xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
       ylabel('Phase interface, {\itx}_f (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
       legend('Analytical','FDS min','FDS max','Location','SouthEast')
    else
       h=plot(t,T(:,1),'k-',t,T(:,2),'r-',t,T(:,3),'g-'); hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,2),'kd');hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,3),'ro');hold on
       h=plot(M_fds(1:5:nt_fds,1),M_fds(1:5:nt_fds,4),'gs');
       xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
       ylabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
       legend('Analytical 1 cm','Analytical 5 cm','Analytical 10 cm','FDS 1 cm','FDS 5 cm','FDS 10 cm')
    end

    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    % Plot attributes

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    % add Git revision if file is available

    Git_Filename = ['pcm_slab_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    % Create the PDF files

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    if i_plot==1
        print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/pcm_slab_xf')
    else
        print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/pcm_slab_T')
    end

end
