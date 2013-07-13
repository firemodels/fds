% McDermott
% 9-22-08
% energy_decay.m
%
% Plot energy decay for Comte-Bellot and Corrsin data
%
% [] = energy_decay(chid,N)
%
% chid = CHID from FDS input file
% N = number of cells in 1D
%
% Example>> energy_decay('csmag_32',32)

function [] = energy_decay(chid,N)

close all
addpath('../../Verification/Turbulence')

L = 0.56549; % \approx 9*(2*pi)/100;

% Gather the Comte-Bellot/Corrsin data
M = load('cbcdata.txt');
k  = M(:,1)*1e2;
E1 = M(:,2)/1e6;
E2 = M(:,3)/1e6;
E3 = M(:,4)/1e6;

% Apply transfer function
k0 = 2*pi/L;
kc = 0.5*N*k0;
delta = pi/kc;
%G = sin(1/2*k.*delta)./(1/2*k.*delta);  % box filter
G = ones(length(k),1);
for j = 1:length(k)
  if k(j)>kc
     tmp = 3*kc/k(j) - 2;
     if tmp > 0
        G(j) = sqrt(tmp); % RJM filter
        %G(j) = 0; % spectral cutoff
     else
        G(j) = 0;
     end
  end
end

E1 = G.*G.*E1;
E2 = G.*G.*E2;
E3 = G.*G.*E3;

% % Use this to check that spectra are correct
% figure(2)
% loglog(k,E1);hold on
% loglog(k,E2)
% loglog(k,E3)
% xlabel('\kappa (1/m)')
% ylabel('E(\kappa) (m^3/s^2)')

% Now integrate
E1_bar = 0;
E2_bar = 0;
E3_bar = 0;
for j = 1:length(k)-1
   dk = k(j+1)-k(j);
   %E1_bar = E1_bar + sqrt(E1(j)*E1(j+1))*dk;
   %E2_bar = E2_bar + sqrt(E2(j)*E2(j+1))*dk;
   %E3_bar = E3_bar + sqrt(E3(j)*E3(j+1))*dk;
   E1_bar = E1_bar + 0.5*(E1(j)+E1(j+1))*dk;
   E2_bar = E2_bar + 0.5*(E2(j)+E2(j+1))*dk;
   E3_bar = E3_bar + 0.5*(E3(j)+E3(j+1))*dk;
end

figure(1)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

if strcmp(chid,'csmag_32')
    if ~exist('mu0_32_devc.csv')
        display('Error: File mu0_32_devc.csv does not exist. Skipping case.')
        return
    end
    M = csvread('mu0_32_devc.csv',2,0);
    t = M(:,1);
    ke = M(:,2);
    H(1)=plot(t,ke,'k-.'); hold on
    
    if ~exist('csmag0_32_devc.csv')
        display('Error: File csmag0_32_devc.csv does not exist. Skipping case.')
        return
    end
    M = csvread('csmag0_32_devc.csv',2,0);
    t = M(:,1);
    ke = M(:,2);
    H(2)=plot(t,ke,'k--'); hold on
end

if ~exist([chid,'_devc.csv'])
    display(['Error: File ' [chid,'_devc.csv'] ' does not exist. Skipping case.'])
    return
end
M = csvread([chid,'_devc.csv'],2,0);
t = M(:,1);
ke = M(:,2);
H(3)=plot(t,ke,'k-'); hold on

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Key_Font_Size)

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Kinetic Energy (m^2/s^2)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

H(4)=plot(0.0,E1_bar,'ko','MarkerSize',8);
plot(0.28,E2_bar,'ro','MarkerSize',8)
plot(0.66,E3_bar,'bo','MarkerSize',8)

if strcmp(chid,'csmag_32')
    legend(H(1:4),'FDS zero visc','FDS mol visc','FDS eddy visc','Filtered CBC data','Location','East')
else
    legend(H(3:4),'FDS','Filtered CBC data','Location','NorthEast')
end

% add SVN if file is available

SVN_Filename = [chid,'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',chid,'_decay'])
