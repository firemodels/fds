% McDermott
% 1/21/2022
% blasius_thermal.m
%
% solution by Pohlhausen, 1921.

close all
clear all

datdir = '../../../out/Convection/';
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

chid = {'Pohlhausen_Pr_p5','Pohlhausen_Pr_1','Pohlhausen_Pr_2'};

% blasius profile

Pr=[0.5, 1, 2]; % Prandtl numbers
u0=1;        % see FDS input file
zmax=1;      % see FDS input file
mu = 0.001;  % see FDS input and output files
rho = 1.2;   % taken from FDS output file
nu = mu/rho;
x = 5; % corresponds to the position of the DEVC array in the FDS input file
neta=500;
[eta,f,fp] = blasius_analytic(u0, zmax, mu, rho, x, neta, .3318);
etamax=eta(end);

% thermal profile from Pohlhausen solution, Holman eq. (B-14)

f1_key={};
f2_key={};

plot_style
Font_Interpreter='LaTeX';
figure

for p=1:length(Pr)

    deta = eta(2)-eta(1);

    I = find(eta<=etamax);
    for i=1:length(I)
        ieta=I(i);
        num=0;
        for j=1:ieta
            intf=trapz(f(1:j))*deta;
            num=num+exp(-Pr(p)/2*intf)*deta;
        end
        denom=num;
        for j=ieta+1:neta
            intf=trapz(f(1:j))*deta;
            denom=denom+exp(-Pr(p)/2*intf)*deta;
        end
        theta(ieta)=(num-deta)/(denom-deta);
    end

    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    H(p)=plot(eta(I),1-theta); hold on
    f1_key{end+1}=['Pr=',num2str(Pr(p))];
    xlabel('$\eta=z/\sqrt{\nu x/u_0}$','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('$\frac{T-T_0}{T_w-T_0}=1-\theta$','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    axis([0 8 0 1])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

end

set(gcf,'Visible',Figure_Visibility);
lh=legend(H,f1_key);

git_file = [datdir, 'Pohlhausen_Pr_2_git.txt'];
addverstr(gca,git_file,'linear')

set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Pohlhausen_similarity_solution']);

% now plot the temperature profile in physical space

figure

for p=1:length(Pr)

    Tw=21; % wall temerature
    T0=20; % ambient temperature

    T = Tw*ones(1,length(theta))+(T0-Tw)*theta;
    z = eta(I)*sqrt(nu*x/u0);

    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    f2_key{end+1}=['Pohlhausen Pr=',num2str(Pr(p))];
    H2(length(f2_key))=plot(T,z); hold on
    line_color=get(H2(end),'Color');

    % add FDS results

    M = importdata([datdir,char(chid(p)),'_line.csv'],',',2);
    zfds=M.data(:,find(strcmp(M.colheaders,'z')));
    Tfds=M.data(:,find(strcmp(M.colheaders,'Tz')));
    Ufds=M.data(:,find(strcmp(M.colheaders,'Uz')));

    f2_key{end+1}=['FDS Pr=',num2str(Pr(p))];
    H2(length(f2_key))=plot(Tfds(2:end),zfds(2:end),'--*','Color',line_color); hold on

    xlabel('$T$ ($^\circ$C)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('$z$ (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    axis([T0-.1 Tw+.1 0 zmax])

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

end

set(gcf,'Visible',Figure_Visibility);
lh2=legend(H2,f2_key);

git_file = [datdir, 'Pohlhausen_Pr_2_git.txt'];
addverstr(gca,git_file,'linear')

set(lh2,'FontName',Font_Name,'FontSize',Key_Font_Size)
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'Pohlhausen_Tz_profile']);

% % plot velocity field (debug)

% figure(3)

% z_bla=eta*sqrt(nu*x/u0);
% u_bla=fp*u0;

% plot(u_bla,z_bla); hold on
% plot(Ufds,zfds)
% axis([0 1.4*u0 0 zmax])
% xlabel('u (m/s)')
% ylabel('z (m)')
% title('Blasius velocity profile')







