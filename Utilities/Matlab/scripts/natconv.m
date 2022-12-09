% McDermott and Henderson
% 26 Oct 2017
% natconv.m

close all
clear all

plot_style

results_dir = ['../../../out/Convection/'];

% calculations below were used for input file setup

g = 9.8;
T1 = [295   303   333   503];
T2 = 293;
Tm = 0.5*(T1+T2);
beta = 1./Tm;
MW = 28.85476; % FDS 'LJ AIR'
P0 = 101325;
rho = P0*MW./(8341.5*Tm);
mu = 1.8216e-5;
cp = 1000;
k=0.018216; % for Pr=1 fluid

Pr = cp*mu/k;
nu = mu./rho;
alpha = k./(rho*cp);

setup=1;

if setup
    % linear region

    % see J.P. Holman p. 361 for correlations

    S0 = 0.002;
    H0 = 16*S0;

    Ra0 = (g*beta.*(T1-T2)*S0^3)./(alpha.*nu);

    Nu0 = ones(1,length(Ra0));
    Tau0 = 1./Nu0 * S0^2./alpha;

    % laminar boundary layer

    S1 = 0.02;
    H1 = 16*S1;

    Ra1 = (g*beta.*(T1-T2)*S1^3)./(alpha.*nu);

    Nu1 = 0.197*Ra1.^(.25)*(H1/S1)^(-1/9);
    Tau1 = 1./Nu1 * S1^2./alpha;

    % turbulent boundary layer

    S2 = .2;
    H2 = 16*S2;

    Ra2 = (g*beta.*(T1-T2)*S2^3)./(alpha.*nu);

    Nu2 = 0.073*Ra2.^(1/3)*(H2/S2)^(-1/9);
    Tau2 = 1./Nu2 * S2^2./alpha;

    S3 = 2;
    H3 = 16*S3;

    Ra3 = (g*beta.*(T1-T2)*S3^3)./(alpha.*nu);

    Nu3 = 0.073*Ra3.^(1/3)*(H3/S3)^(-1/9);
    Tau3 = 1./Nu3 * S3^2./alpha;

    S4 = 20;
    H4 = 16*S4;

    Ra4 = (g*beta.*(T1-T2)*S4^3)./(alpha.*nu);

    Nu4 = 0.073*Ra4.^(1/3)*(H4/S4)^(-1/9);
    Tau4 = 1./Nu4 * S4^2./alpha;

    % % high Ra (Nature paper)

    % Ra4 = logspace(7,17,11);
    % Nu4 = 0.124*Ra3.^(0.309)*(H2/S2)^(-1/9);

    Ra_Limit_1 = 2000;
    Ra_Limit_2 = 6000;
    Ra_Limit_3 = 2e5;

    RAYLEIGH = logspace(0,14,1000);
    for i=1:length(RAYLEIGH)
        if RAYLEIGH(i)<Ra_Limit_1
            NUSSELT(i)=1;
        elseif RAYLEIGH(i)>Ra_Limit_1 & RAYLEIGH(i)<=Ra_Limit_2
            NUSSELT(i)=0.197*RAYLEIGH(i)^(.25) * (16)^(-1/9);
        elseif RAYLEIGH(i)>Ra_Limit_2 & RAYLEIGH(i)<=Ra_Limit_3
            NUSSELT(i)=0.197*RAYLEIGH(i)^(.25) * (16)^(-1/9);
        elseif RAYLEIGH(i)>Ra_Limit_3
            NUSSELT(i)=0.073*RAYLEIGH(i)^(1/3) * (16)^(-1/9) ;
        end
    end

    figure(1)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    marker_handle(1)=loglog(RAYLEIGH,NUSSELT,'k-'); hold on
    % loglog(Ra1,Nu1,'k-')
    % loglog(Ra2,Nu2,'k-')
    % loglog(Ra3,Nu3,'k-')
    % loglog(Ra4,Nu4,'k-')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    axis([1 1e15 .5 1e4])
    xlabel('Rayleigh Number','FontSize',Label_Font_Size)
    ylabel('Nusselt Number','FontSize',Label_Font_Size)
    %title('Natural Convection in Vertical Enclosures')

    %axis([min(Ra0),max(Ra3),.5,max(Nu3)])

    % th=annotation('textarrow',[.35,.35],[.4,.25],'String','Nu_s = 0.197 Ra_s^{1/4} (H/s)^{-1/9}');
    % set(th,'FontSize',16)

    % th=annotation('textarrow',[.542,.63],[.55,.55],'String','Nu_s = 0.073 Ra_s^{1/3} (H/s)^{-1/9}');
    % set(th,'FontSize',16)

    % return
end

% FDS results

casename={...
'natconv_1',...
'natconv_2',...
'natconv_3',...
'natconv_4',...
'natconv_5',...
'natconv_6',...
'natconv_7',...
'natconv_8',...
'natconv_9',...
'natconv_10',...
'natconv_11',...
'natconv_12',...
'natconv_13',...
'natconv_14',...
'natconv_15',...
'natconv_16',...
'natconv_17'...
};

delta = [0.002 0.02 0.02 0.02 0.02 0.2 0.2 0.2 0.2 2 2 2 2 20 20 20 20];
T = [303 295 303 333 503 295 303 333 503 295 303 333 503 295 303 333 503];

marker_style = {'r^','bsq','k+'};
res = {'8','16','32'};
qcolhdrs = {'Q1-1','Q1-2','Q1-3','Q1-4','Q1-5','Q2-1','Q2-2','Q2-3','Q2-4','Q2-5'};

for j=1:length(res)
    for i=1:length(delta)

        M = importdata([results_dir,casename{i},'_',res{j},'_devc.csv']);

        % check for steady state
        t = M.data(:,1);
        Q1 = M.data(:,2);
        Q2 = M.data(:,3);
        rho = M.data(end,5);
        alpha = k/(rho*cp);
        nu = mu/rho;
        b = 2./(T(i)+T2);
        Ra(i) = (g*b*(T(i)-T2)*delta(i)^3)/(alpha*nu);

        % figure(2)
        % plot(t,Q1,'r-'); hold on
        % plot(t,Q2,'b-');
        % plot(t,Q1+Q2,'k--')

        M = importdata([results_dir,casename{i},'_',res{j},'_line.csv'],',',2);

        col = [];
        for n=1:length(qcolhdrs)
            col = [col,find(strcmp(M.colheaders,qcolhdrs{n}))];
        end
        col;

        Q = mean(mean(abs(M.data(:,col))*1000));  % heat flow, W

        Nu_FDS = Q*(delta(i)/k)/(T(i)-T2);

        % figure(1)
        set(gca,'Units',Plot_Units)
        set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
        marker_handle(j+1)=plot(Ra(i),Nu_FDS,marker_style{j},'MarkerSize',8);
    end
end

lh=legend(marker_handle,'Correlation, Nu = C Ra^n (L/H)^m','FDS {\itH/\Deltax}=8','FDS {\itH/\Deltax}=16','FDS {\itH/\Deltax}=32','Location','Northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'natconv_17_32_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natural_convection_vertical_enclosure');


% Complex geometry cases (rotated grid 18 degrees)

if 1

setup=1;

if setup
    % linear region

    % see J.P. Holman p. 361 for correlations

    S0 = 0.002;
    H0 = 16*S0;

    Ra0 = (g*beta.*(T1-T2)*S0^3)./(alpha.*nu);

    Nu0 = ones(1,length(Ra0));
    Tau0 = 1./Nu0 * S0^2./alpha;

    % laminar boundary layer

    S1 = 0.02;
    H1 = 16*S1;

    Ra1 = (g*beta.*(T1-T2)*S1^3)./(alpha.*nu);

    Nu1 = 0.197*Ra1.^(.25)*(H1/S1)^(-1/9);
    Tau1 = 1./Nu1 * S1^2./alpha;

    % turbulent boundary layer

    S2 = .2;
    H2 = 16*S2;

    Ra2 = (g*beta.*(T1-T2)*S2^3)./(alpha.*nu);

    Nu2 = 0.073*Ra2.^(1/3)*(H2/S2)^(-1/9);
    Tau2 = 1./Nu2 * S2^2./alpha;

    S3 = 2;
    H3 = 16*S3;

    Ra3 = (g*beta.*(T1-T2)*S3^3)./(alpha.*nu);

    Nu3 = 0.073*Ra3.^(1/3)*(H3/S3)^(-1/9);
    Tau3 = 1./Nu3 * S3^2./alpha;

    S4 = 20;
    H4 = 16*S4;

    Ra4 = (g*beta.*(T1-T2)*S4^3)./(alpha.*nu);

    Nu4 = 0.073*Ra4.^(1/3)*(H4/S4)^(-1/9);
    Tau4 = 1./Nu4 * S4^2./alpha;

    % % high Ra (Nature paper)

    % Ra4 = logspace(7,17,11);
    % Nu4 = 0.124*Ra3.^(0.309)*(H2/S2)^(-1/9);

    Ra_Limit_1 = 2000;
    Ra_Limit_2 = 6000;
    Ra_Limit_3 = 2e5;

    RAYLEIGH = logspace(0,14,1000);
    for i=1:length(RAYLEIGH)
        if RAYLEIGH(i)<Ra_Limit_1
            NUSSELT(i)=1;
        elseif RAYLEIGH(i)>Ra_Limit_1 & RAYLEIGH(i)<=Ra_Limit_2
            NUSSELT(i)=0.197*RAYLEIGH(i)^(.25) * (16)^(-1/9);
        elseif RAYLEIGH(i)>Ra_Limit_2 & RAYLEIGH(i)<=Ra_Limit_3
            NUSSELT(i)=0.197*RAYLEIGH(i)^(.25) * (16)^(-1/9);
        elseif RAYLEIGH(i)>Ra_Limit_3
            NUSSELT(i)=0.073*RAYLEIGH(i)^(1/3) * (16)^(-1/9) ;
        end
    end

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    marker_handle(1)=loglog(RAYLEIGH,NUSSELT,'k-'); hold on
    % loglog(Ra1,Nu1,'k-')
    % loglog(Ra2,Nu2,'k-')
    % loglog(Ra3,Nu3,'k-')
    % loglog(Ra4,Nu4,'k-')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    axis([1 1e15 .5 1e4])
    xlabel('Rayleigh Number','FontSize',Label_Font_Size)
    ylabel('Nusselt Number','FontSize',Label_Font_Size)
    %title('Natural Convection in Vertical Enclosures')

    %axis([min(Ra0),max(Ra3),.5,max(Nu3)])

    % th=annotation('textarrow',[.35,.35],[.4,.25],'String','Nu_s = 0.197 Ra_s^{1/4} (H/s)^{-1/9}');
    % set(th,'FontSize',16)

    % th=annotation('textarrow',[.542,.63],[.55,.55],'String','Nu_s = 0.073 Ra_s^{1/3} (H/s)^{-1/9}');
    % set(th,'FontSize',16)

    % return
end

% FDS results

casename={...
'natconv_1',...
'natconv_2',...
'natconv_3',...
'natconv_4',...
'natconv_5',...
'natconv_6',...
'natconv_7',...
'natconv_8',...
'natconv_9',...
'natconv_10',...
'natconv_11',...
'natconv_12',...
'natconv_13',...
'natconv_14',...
'natconv_15',...
'natconv_16',...
'natconv_17'...
};

delta = [0.002 0.02 0.02 0.02 0.02 0.2 0.2 0.2 0.2 2 2 2 2 20 20 20 20];
T = [303 295 303 333 503 295 303 333 503 295 303 333 503 295 303 333 503];

marker_style = {'r^','bsq'};
res = {'8','16'};

for j=1:length(res)
    for i=1:length(delta)

        M = importdata([results_dir,casename{i},'_',res{j},'_rot_18_devc.csv']);

        % check for steady state
        t = M.data(:,1);
        Q1 = M.data(:,find(strcmp(M.colheaders,'"Q1-0"')));
        Q2 = M.data(:,find(strcmp(M.colheaders,'"Q2-0"')));
        rho = M.data(end,find(strcmp(M.colheaders,'"rho"')));
        alpha = k/(rho*cp);
        nu = mu/rho;
        b = 2./(T(i)+T2);
        Ra(i) = (g*b*(T(i)-T2)*delta(i)^3)/(alpha*nu);

        % figure(2)
        % plot(t,Q1,'r-'); hold on
        % plot(t,Q2,'b-');
        % plot(t,Q1+Q2,'k--')

        qrange = find(t>t(end)/2);
        A = (delta(i)*4)*(delta(i)*16);
        q = mean(abs(Q2(qrange)))*1000/A;  % heat flux, W/m2

        Nu_FDS = q*(delta(i)/k)/(T(i)-T2);

        % figure(1)
        set(gca,'Units',Plot_Units)
        set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
        marker_handle(j+1)=plot(Ra(i),Nu_FDS,marker_style{j},'MarkerSize',8);
    end
end

lh=legend(marker_handle(1:3),'Correlation, Nu = C Ra^n (L/H)^m','GEOM {\itH/\Deltax}=8','GEOM {\itH/\Deltax}=16','Location','Northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'natconv_17_16_rot_18_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natconv_geom');

end
