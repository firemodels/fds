% McDermott
% 2-23-2021
% spray_pattern.m
%
% help clarify spray_pattern parameters
% plots are saved to fds/Manuals/FDS_User_Guide/FIGURES/spray_pattern*

close all
clear all
plot_style

b = [0, 5, 10, 100, 1000];
mu = [0,10,20]  * pi/180;
mu_str = {'0' '10' '20'};

phi_min = 0;
phi_max = [22.5, 45, 90] * pi/180;
phi_str = {'22p5' '45' '90'};

for k=1:length(phi_max)
    phi = linspace(phi_min,phi_max(k),1000);
    for j=1:length(mu)
        figure
        set(gca,'Units',Plot_Units)
        set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
        for i=1:length(b)
            f = exp(-b(i)*((phi-mu(j))/(phi_max(k)-phi_min)).^2);
            H(i)=plot(phi*180./pi,f); hold on
            Key{i}=['\beta=',num2str(b(i))];
        end
        axis([0, 90, 0, 1.1]);
        set(gca, 'YDir','reverse');
        xlabel('\phi (degrees)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size);
        ylabel('{\itf} (\phi)'       ,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size);
        set(gca,'FontName',Font_Name);
        set(gca,'FontSize',Label_Font_Size);
        lh=legend(H,Key,'location','southeast');
        set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size);
        th=text(45,0.1,['\mu=',num2str(mu(j)*180/pi),', \phi_{max}=',num2str(phi_max(k)*180/pi),' degrees']);
        set(th,'FontName',Font_Name,'FontSize',Label_Font_Size);
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'Units',Paper_Units);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
        print(gcf,'-dpdf',['../../../Manuals/FDS_User_Guide/FIGURES/spray_pattern_mu_',mu_str{j},'_phimax_',phi_str{k}]);
        clear H Key;
    end
end





