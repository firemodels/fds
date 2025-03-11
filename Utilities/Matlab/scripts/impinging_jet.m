% McDermott
% 12-8-2023
% impinging_jet.m

close all
clear all

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

outdir = '../../../out/Convection/';
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
chid = 'impinging_jet';

% plot the correlation

D_h = 0.2;                            % hydraulic diameter of the jet [m]
H   = 1;                              % distance from jet to wall [m]
mu  = 1.822E-05;                      % dynamic viscosity, [kg/m/s]
k   = 2.566E-02;                      % thermal conductivity [W/m/K]
Pr  = 0.71;                           % Prandtl number
T_j = 100;                            % tempetature of jet exit [C]
T_w = 20;                             % constant plate temperature [C]
rho_a = 1.200;                        % ambient density [kg/m3]
rho_j = rho_a * (T_w+273)/(T_j+273);  % jet fluid density [kg/m3]
nu = mu/rho_j;                        % kinematic viscosity [m2/s]
U_j = [10,40];                        % jet exit velocity [m/s]
A_r  = 0.01;                          % D_h^2/(4 r^2), where r defines the outer extent of the averaging region

% Plot correlation versus Re_j for a given geometry

% Martin corrleation
G = 2*sqrt(A_r)*( (1-2.2*sqrt(A_r))/(1+0.2*(H/D_h-6)*sqrt(A_r)) );

RE = linspace(2e3,5e5);
NU = G*2*sqrt(RE).*sqrt(1+0.005*RE.^0.55).* Pr^0.42;

K(1)=plot(RE,NU,'k-','linewidth',2); hold on

res_str = {'coarse','medium','fine'};
Re_str  = {'1e5','4e5'};

% relative error tolerance
E_tol = 0.1;

E_FDS = zeros(3,2);

for j=1:length(res_str)
    for i=1:length(Re_str)

        Re_j = D_h*U_j(i)/nu;  % jet Reynolds number

        Nu = G*2*sqrt(Re_j)*sqrt(1+0.005*Re_j^0.55) * Pr^0.42;

        h_cor = Nu * k/D_h;

        % FDS results

        M = importdata([outdir,'impinging_jet_Re_',Re_str{i},'_',res_str{j},'_devc.csv'],',',2);
        HF = mean(M.data(floor(end/2):end,find(strcmp(M.colheaders,'"HF"'))));

        qconv = HF * 1000; % W/m2
        h_fds = qconv/(T_j-T_w);
        Nu_fds = h_fds*D_h/k;

        E_FDS(i,j) = abs(Nu - Nu_fds)/abs(Nu);

        if E_FDS(i,j) > E_tol
          disp(['Matlab Warning: impinging jet error = ',num2str(E_FDS(i,j)),' at Re_j=',Re_str{i},', Res=',res_str{j}])
        end

        if i==1
            if j==1
                K(2)=plot(Re_j,Nu_fds,'bsq','linewidth',2); Key={'FDS coarse'};
            elseif j==2
                K(3)=plot(Re_j,Nu_fds,'rsq','linewidth',2); Key={'FDS coarse','FDS medium'};
            elseif j==3
                K(4)=plot(Re_j,Nu_fds,'gsq','linewidth',2); Key={'FDS coarse','FDS medium','FDS fine'};
            end
        elseif i==2
            if j==1
                plot(Re_j,Nu_fds,'bsq','linewidth',2)
            elseif j==2
                plot(Re_j,Nu_fds,'rsq','linewidth',2)
            elseif j==3
                plot(Re_j,Nu_fds,'gsq','linewidth',2)
            end
        end

    end
end

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Re','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Nu','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend(K,['Martin',Key],'location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size);

Git_Filename = [outdir,'impinging_jet_Re_1e5_coarse_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'impinging_jet_correlation']);


% plot the profile of heat transfer coefficient

style = {'k-.','k--','k-'};
clear K lh;

for i=1:length(Re_str)
    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    for j=1:length(res_str)

        M = importdata([outdir,'impinging_jet_Re_',Re_str{i},'_',res_str{j},'_line.csv'],',',2);
        x = M.data(:,1);
        q_x = M.data(:,find(strcmp(M.colheaders,'QCONV'))) * 1000;
        Nu_x = q_x/(T_j-T_w)*D_h/k;

        K(j)=plot(x,Nu_x,style{j}); hold on
    end
    xlabel('{\itx} (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('Nu_{Dh}({\itx})','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

    axis([-.5 .5 0 1000])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    lh=legend(K,'D_h/\delta{\itx}=7','D_h/\delta{\itx}=14','D_h/\delta{\itx}=28','location','northwest');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size);

    Git_Filename = [outdir,'impinging_jet_Re_1e5_coarse_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,'impinging_jet_local_Re_',Re_str{i}]);

    hold off
end





