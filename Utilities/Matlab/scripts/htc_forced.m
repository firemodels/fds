% McDermott
% 11-14-2022
% htc_forced.m

close all
clear all

plot_style

outdir = '../../../out/Convection/';
pltdir = '../../../fds/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
chid = 'forced_conv_flat_plate';
vel = {'u0p1','u1','u10'};
res = {'25cm','10cm','2p5cm'};
res_style = {'b*','ksq','mo'};
u = [0.1,1,10];
dx = [0.25,0.1,0.025]; % grid resolution (m)
dx_style = {'b--','k--','m--'};

N = 32; % number of points in line device
L = 16; % domain length
rho = 1.165; % pure N2 specified in FDS input file
mu = 1E-05;  % specified in FDS input file
k = 1E-02;   % specified in FDS input file
T_w = 30;    % constant wall temperature specified
T_g = 20;    % T_infty inlet temperature specified

x = linspace(0,L,N);

for i=1:length(vel)

    Re_x = (rho/mu)*u(i)*x;
    Re_1 = (rho/mu)*u(i)*1;

    Nu_x = 0.0296 * Re_x.^0.8; % Incropera and Dewitt correlation, Eq. 7.36, Table 7.7
    Nu_1 = 0.0296 * Re_1.^0.8; % reference value for normalization of errors

    % BL_thickness = 0.37*x.*Re_x.^(-0.2);
    % BL_thickness(end)

    figure(i)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    H(1)=plot(Re_x,Nu_x); hold on
    xlabel('Re_{\itx}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('Nu_{\itx}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    for j=1:length(res)

        % FDS results

        M = importdata([outdir,chid,'_',vel{i},'_',res{j},'_line.csv'],',',2);
        x_fds = M.data(:,find(strcmp(M.colheaders,'HTC-x')));
        h_fds = M.data(:,find(strcmp(M.colheaders,'HTC')));
        q_fds = M.data(:,find(strcmp(M.colheaders,'QCONV')));
        DTMP_fds = -1000*q_fds./h_fds;
        Nu_x_fds = h_fds.*DTMP_fds/(T_w-T_g).*x_fds/k;
        Re_x_fds = (rho/mu)*u(i)*x_fds;

        H(1+j)=plot(Re_x_fds,Nu_x_fds,res_style{j});
        H(1+length(res)+j)=plot(Re_x_fds,2/dx(j)*x_fds,dx_style{j});

        REL_ERROR(i,j) = abs(Nu_x(end)-Nu_x_fds(end))/Nu_x(end);

    end

    xl = get(gca,'XLim');
    yl = get(gca,'YLim');

    xtxt = xl(1) + 0.5*(xl(2)-xl(1));
    ytxt = yl(1) + 0.9*(yl(2)-yl(1));
    text(xtxt,ytxt,['Velocity = ',num2str(u(i)),' m/s'],'FontName',Font_Name,'FontSize',Title_Font_Size)

    lh=legend(H,'forced local','FDS {\it\deltax}=25cm','FDS {\it\deltax}=10cm','FDS {\it\deltax}=2.5cm','2*1/0.25*{\itx}','2*1/0.10*{\itx}','2*1/0.025*{\itx}','location','northwest');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    Git_Filename = [outdir,chid,'_',vel{i},'_25cm_git.txt'];
    addverstr(gca,Git_Filename,'linear')

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,chid,'_',vel{i}]);
end

% Error test

% REF_ERROR determined from initial run
%            dx=25cm   dx=10cm   dx=2.5cm
REF_ERROR = [0.1010    0.1930    0.0235; ...
             0.2628    0.0330    0.2927; ...
             0.3371    0.1494    0.1421];

% Conclusion: Forced convection has approx 30% error.

ERROR = norm(REL_ERROR-REF_ERROR); % 0.1 tolerance

if ERROR > 0.1
   display(['Matlab Warning: Forced convection out of tolerance. ERROR = ',num2str(ERROR)])
end
