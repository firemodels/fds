% McDermott
% 19 Sep 2017
% McCaffrey_Correlation.m

close all
clear all

expdir = '../../../exp/Submodules/macfp-db/Gaseous_Pool_Fires/McCaffrey_Flames/Computational_Results/2017/Data/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/McCaffrey_Plume/';

plot_style

figure(1)
reset(gcf)
reset(gca)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

figure(2)
reset(gcf)
reset(gca)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])


Q = [14.4 21.7 33.0 44.9 57.5]; % kW
g = 9.8; % m/s2
rho = 1.18; % kg/m3
cp = 1; % kJ/(kg*K)
T0 = 273.15 + 20; % K
D = 0.3; % m

QS = (Q/(rho*cp*T0*sqrt(g)*D^(5/2)));
DS = (Q/(rho*cp*T0*sqrt(g))).^(2/5); % m

chid = {'McCaffrey_14kW','McCaffrey_22kW','McCaffrey_33kW','McCaffrey_45kW','McCaffrey_57kW'};
mark = {'ko','k+','k^','ksq','kd'};
n_chid = length(chid);

% McCaffrey plume correlations NBSIR 79-1910

zq = logspace(-2,0,100);

for i=1:length(zq)
    if zq(i)<0.08
        vq(i) = 6.84*zq(i)^0.5;
        Tq(i) = 800*zq(i)^0;
    elseif zq(i)>=0.08 & zq(i)<=0.2
        vq(i) = 1.93*zq(i)^0;
        Tq(i) = 63*zq(i)^(-1);
    elseif zq(i)>0.2
        vq(i) = 1.12*zq(i)^(-1/3);
        Tq(i) = 21.6*zq(i)^(-5/3);
    end
end

% % write correlations to a file

% % centerline temperature

% fid = fopen('../Experimental_Data/McCaffrey_Correlation_Temperature.csv','wt');
% fprintf(fid,'%s, %s\n','Z/Q^0.4','dT (C)');
% for i=1:length(zq)
%     fprintf(fid,'%f, %f\n',zq(i),Tq(i));
% end
% fclose(fid);

% % centerline velocity

% fid = fopen('../Experimental_Data/McCaffrey_Correlation_Velocity.csv','wt');
% fprintf(fid,'%s, %s\n','Z/Q^0.4','V/Q^0.2');
% for i=1:length(zq)
%     fprintf(fid,'%f, %f\n',zq(i),vq(i));
% end
% fclose(fid);

h1 = zeros(1,length(chid)+1);
h2 = zeros(1,length(chid)+1);

figure(1); h1(end)=loglog(zq,vq,'b--','linewidth',2); hold on
figure(2); h2(end)=loglog(zq,Tq,'r--','linewidth',2); hold on

for i=1:length(chid)
    V = importdata([expdir,chid{i},'_V.csv'],',',1);
    T = importdata([expdir,chid{i},'_T.csv'],',',1);
    figure(1); h1(i)=loglog(V.data(:,1),V.data(:,2),mark{i});
    figure(2); h2(i)=loglog(T.data(:,1),T.data(:,2),mark{i});
end

% format and print velocity correlation

figure(1)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('{\itz/Q}^{2/5}','FontSize',Label_Font_Size)
ylabel('{\itV/Q}^{1/5}','FontSize',Label_Font_Size)

xmin = 0.01;
xmax = 1;
ymin = .5;
ymax = 3.0;
axis([xmin xmax ymin ymax])

lh = legend(h1,'14.4 kW','21.7 kW','33.0 kW','44.9 kW','57.5 kW','({\itz/Q}^{2/5})^\eta','Location','SouthEast');
set(lh,'FontSize',Key_Font_Size)

text(.0125,2.6,'McCaffrey Centerline Velocity Correlation','FontName',Font_Name,'FontSize',Title_Font_Size)
text(.04,1.2,'\eta=1/2','FontSize',Label_Font_Size,'FontName',Font_Name)
text(.10,1.7,'\eta=0','FontSize',Label_Font_Size,'FontName',Font_Name)
text(.295,1.25,'\eta=-1/3','FontSize',Label_Font_Size,'FontName',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'McCaffrey_Velocity_Correlation'])

% format and print temperature correlation

figure(2)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('{\itz/Q}^{2/5}','FontSize',Label_Font_Size)
ylabel('\Delta{\itT} (\circC)','FontSize',Label_Font_Size)

axis([.008 1 100 1200])

lh = legend(h2,'14.4 kW','21.7 kW','33.0 kW','44.9 kW','57.5 kW','({\itz/Q}^{2/5})^\eta','Location','SouthWest');
set(lh,'FontSize',Key_Font_Size)

text(.01,1000,'McCaffrey Centerline Temperature Correlation','FontName',Font_Name,'FontSize',Title_Font_Size)
text(.03,650,'\eta=0','FontSize',Label_Font_Size,'FontName',Font_Name)
text(.08,425,'\eta=-1','FontSize',Label_Font_Size,'FontName',Font_Name)
text(.18,125,'\eta=-5/3','FontSize',Label_Font_Size,'FontName',Font_Name)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'McCaffrey_Temperature_Correlation'])

return



% % Baum and McCaffrey plume correlations (in terms of D*)
% % 2nd IAFSS, pp. 129-148

% zs = logspace(-2,1,100);
% n = [1/2 0 -1/3];
% A = [2.18 2.45 3.64];
% B = [2.91 3.81 8.41];

% for i=1:length(zs)
%     if zs(i)<1.32
%         j=1;
%     elseif zs(i)>=1.32 & zs(i)<=3.3
%         j=2;
%     elseif zs(i)>3.3
%         j=3;
%     end
%     us(i) = A(j)*zs(i)^n(j);
%     Ts(i) = B(j)*zs(i)^(2*n(j)-1);
% end
