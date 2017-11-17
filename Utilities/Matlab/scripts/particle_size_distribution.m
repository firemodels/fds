% McDermott
% 23 Oct 2017
% particle_size_distribution.m
%
% This script generates the Rosin-Rammler / log-normal particle size distribution
% plot in the FDS Tech Guide (presently Fig. 8.1).

close all
clear all

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

D_v50 = 1;

gamma = 2.4;
sigma = 1.15/gamma;

N = 100;
D = linspace(0,3,N);
dD = D(2)-D(1);

for i=1:N

    if i==1
        F1(i) = 0;
    else
        DP = 0.5*(D(i-1)+D(i));
        F1(i) = F1(i-1) + exp(-log(DP/D_v50)^2/(2*sigma^2)) * dD / (sqrt(2*pi)*sigma*DP);
    end

    F2(i) = 1 - exp(-0.693*(D(i)/D_v50)^gamma);

    if D(i)<D_v50
        Fv(i) = F1(i);
    else
        Fv(i) = F2(i);
    end

    % probability density function (derivative of CDF)
    if i==1
        fv(i) = 0;
    else
        fv(i) = (Fv(i)-Fv(i-1))/(D(i)-D(i-1));
    end

    if i==1
        Fn(i) = 0;
    else
        DP = 0.5*(D(i-1)+D(i));
        Fn(i) = Fn(i-1) + fv(i)/DP^3 * dD;
    end

end

Fn = Fn/Fn(end);

H(1)=plot(D,F1,'k--');
hold on
H(2)=plot(D,F2,'k-.');
H(3)=plot(D,Fv,'k-','linewidth',2);
H(4)=plot(D,Fn,'b-','linewidth',2);
H(5)=plot([D_v50,D_v50],[0,1],'k:');

xlabel('droplet diameter (mm)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'F_v log-normal','F_v Rosin-Rammler','F_v combined (FDS)','F_n (cumulative number fraction)','D_{v,0.5}',...
    'location','southeast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/particle_size_distribution');








