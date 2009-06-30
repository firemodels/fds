% McDermott
% 3-4-09
% law_of_the_wall.m

close all
clear all

paper_width  = 6.0; % inches
paper_height = 4.5; % inches

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

% figure
% plot(up,zp)
% xlabel('u+')
% ylabel('z+')

figure
font_size = 14;
semilogx(zp,up,'Linewidth',1.5); hold on
semilogy(zp,uu,'r--','Linewidth',1.5)
set(gca,'Units','inches')
set(gca,'FontName','Times')
set(gca,'FontSize',font_size)
set(gca,'Position',[1,0.75,4.25,3.15])
text(1.5,8,'$u^+ = z^+$','Fontsize',font_size,'FontName','Times','Interpreter','LaTeX')
annotation('textarrow',[.5 .6],[.7 .555],'String','$u^+ = 2.4 {\rm ln} z^+ + 5.2$','Fontsize',font_size,'Interpreter','LaTeX');
annotation('textarrow',[.65 .6],[.38 .535],'String','$u^+ = A(z^+)^B$','Fontsize',font_size,'Interpreter','LaTeX');
H = annotation('line',[.34 .34],[.15 .6]);
set(H,'LineStyle','--')
annotation('textarrow',[.44 .34],[.22 .26],'String',' $z^+ = 11.81$','Fontsize',font_size,'Interpreter','LaTeX');
xlabel('$z^+$','Fontsize',font_size,'Interpreter','LaTeX')
ylabel('$u^+$','Fontsize',font_size,'Interpreter','LaTeX')
set(gcf,'Visible','on');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[paper_width paper_height]);
set(gcf,'PaperPosition',[0 0 paper_width paper_height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_Verification_Guide/FIGURES/lawofthewall')

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

