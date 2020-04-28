% McDermott
% 3-4-09
% law_of_the_wall.m

close all
clear all

% plot Werner and Wengle velocity profile

A = 8.3;
B = 1/7;
n = 100;
zp = logspace(0,4,n);
uu = zp; % log law
for j = 1:n
    if zp(j)<11.81
        up(j) = zp(j);
    else
        up(j) = A*zp(j)^B; % Werner and Wengle (power law)
        uu(j) = (1/0.41)*log(zp(j))+5.2; % log law
    end
end

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

semilogx(zp,up,'b-'); hold on
semilogy(zp,uu,'r--')
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
text(1.5,8,'{\itu}^+ = {\itz}^+','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
text(200,15,'{\itu}^+ = 2.4 ln {\itz}^+ + 5.2','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'Color',[1 0 0])
text(500,30,'{\itu}^+ = A({\itz}^+)^B','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'Color',[0 0 1])
line([11.81 11.81],[0 20],'LineStyle','--','Color',[0 0 0])
text(15,5,'{\itz}^+ = 11.81','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter);
xlabel('{\itz}^+','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('{\itu}^+','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/lawofthewall')

% % check FDS parameters
%
% alpha = (1-B)/2*A^((1+B)/(1-B))
% beta = 1+B
% eta = (1+B)/A
% gamma = 2/(1+B)
% H = 1;
% N = 16;
% dz = H/N
% rho = 1.198;
% mu = 1.81e-5;
% nu = mu/rho;
% nu_over_dz = nu/dz
%
% u1 = 8.311075e-2;
% tau_w = ( alpha*nu_over_dz^beta + eta*nu_over_dz^B*abs(u1) )^gamma
% zplus = dz*sqrt(tau_w)/nu

