% McDermott
% 21 Sep 2017
% McCaffrey_Surface_Temp.m
%
% Extrapolate McCaffrey TC data to burner surface for FDS boundary condition.

close all
clear all

expdir = '../../../exp/Submodules/macfp-db/Gaseous_Pool_Fires/McCaffrey_Flames/Computational_Results/2017/Data/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/McCaffrey_Plume/';

plot_style

Q = [14.4 21.7 33.0 44.9 57.5]; % kW
g = 9.8; % m/s2
rho = 1.18; % kg/m3
cp = 1; % kJ/(kg*K)
T0 = 273.15 + 20; % K
D = 0.3; % m

QS = (Q/(rho*cp*T0*sqrt(g)*D^(5/2)));
DS = (Q/(rho*cp*T0*sqrt(g))).^(2/5); % m

chid = {'McCaffrey_14kW','McCaffrey_22kW','McCaffrey_33kW','McCaffrey_45kW','McCaffrey_57kW'};
key  = {'14.4 kW','21.7 kW','33.0 kW','44.9 kW','57.5 kW'};
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

for i=1:length(chid)

    figure
    reset(gcf)
    reset(gca)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    h1 = zeros(1,3);
    h1(2)=loglog(zq,Tq,'r--','linewidth',2); hold on

    T = importdata([expdir,chid{i},'_T.csv'],',',1);
    ZS = T.data(:,1);
    TC = T.data(:,2);
    h1(1)=loglog(ZS,TC,mark{i});

    % identify and color data near surface
    isurf = find(ZS<0.05);
    loglog(ZS(isurf),TC(isurf),'bo','MarkerSize',10,'MarkerFaceColor','b')

    % least squares fit to near surface data
    A = [ZS(isurf),ones(length(isurf),1)];
    b = TC(isurf);
    x = pinv(A)*b; % solves linear least squares problem
    T_SURF = x(2);

    % plot best fit line
    z1 = 0.008;
    z2 = 0.05;
    T1 = x(1)*z1 + x(2);
    T2 = x(1)*z2 + x(2);
    h1(3)=loglog([z1 z2],[T1 T2],'g-','LineWidth',3);

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    xlabel('{\itz/Q}^{2/5}','FontSize',Label_Font_Size)
    ylabel('\Delta{\itT} (\circC)','FontSize',Label_Font_Size)

    axis([.008 1 100 1200])

    lh = legend(h1,key{i},'({\itz/Q}^{2/5})^\eta','Least Squares Fit','Location','SouthWest');
    set(lh,'FontSize',Key_Font_Size)

    text(.01,1000,['McCaffrey Centerline Temperature Data, ',key{i}],'FontName',Font_Name,'FontSize',Title_Font_Size)
    text(.03,650,'\eta=0','FontSize',Label_Font_Size,'FontName',Font_Name)
    text(.08,425,'\eta=-1','FontSize',Label_Font_Size,'FontName',Font_Name)
    text(.18,125,'\eta=-5/3','FontSize',Label_Font_Size,'FontName',Font_Name)

    text(.009,400,['{\itT}_{SURF}=',num2str(T_SURF,'%4.0f'),' \circC'],'FontSize',Label_Font_Size,'FontName',Font_Name)

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,chid{i},'_Surface_Temp'])

    clf
    hold off
end

