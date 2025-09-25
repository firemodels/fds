% McDermott and Moon
% 8-15-12
% heated_channel.m
%
% References:
%
% Moser, Kim, and Mansour. DNS of Turbulent Channel Flow up to Re_tau=590.
% Physics of Fluids, vol 11, 943-945, 1999.
%
% Kim, Moin & Moser. [Numerical Method, Re_tau = 178.12]. J. Fluid Mech.
% vol 177, 133-166, 1987.

close all
clear all

plot_style

% plot the DNS results

outdir='../../../out/Heated_Channel_Flow/';
expdir='../../../exp/Heated_Channel_Flow/';

M = importdata([expdir,'heated_channel_dns_data.csv'],',',1);

yp_up_mean      = M.data(:,find(strcmp(M.colheaders,'yp_up_mean')));
up_mean         = M.data(:,find(strcmp(M.colheaders,'up_mean')));
yp_Tp_mean      = M.data(:,find(strcmp(M.colheaders,'yp_Tp_mean')));
Tp_mean_Pr0p10  = M.data(:,find(strcmp(M.colheaders,'Tp_mean_Pr0p1')));
Tp_mean_Pr0p71  = M.data(:,find(strcmp(M.colheaders,'Tp_mean_Pr0p71')));
Tp_mean_Pr2p00  = M.data(:,find(strcmp(M.colheaders,'Tp_mean_Pr2p0')));
range = 1:3:length(yp_up_mean);

f1=figure;
a1=gca;
set(f1,'Visible',Figure_Visibility);
set(a1,'Units',Plot_Units)
set(a1,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hfig1(2)=semilogx(yp_up_mean,up_mean,'ro'); hold on
axis([1e0 1e3 0 30])

set(a1,'FontName',Font_Name)
set(a1,'FontSize',Title_Font_Size)

xlabel('{\it z}^+','Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('{\it u}^+','Interpreter',Font_Interpreter,'Fontname',Font_Name)


f2=figure;
a2=gca;
set(f2,'Visible',Figure_Visibility);
set(a2,'Units',Plot_Units)
set(a2,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hfig2(2)=semilogx(yp_Tp_mean(range),Tp_mean_Pr0p10(range),'bo'); hold on
semilogx(yp_Tp_mean(range),Tp_mean_Pr0p71(range),'bo')
semilogx(yp_Tp_mean(range),Tp_mean_Pr2p00(range),'bo')
axis([1e0 1e3 0 30])

set(a2,'FontName',Font_Name)
set(a2,'FontSize',Label_Font_Size)

xlabel('{\it z}^+','Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('{\it T}^+','Interpreter',Font_Interpreter,'Fontname',Font_Name)

text(225,26,'Pr=2.0','Interpreter',Font_Interpreter,'Fontname',Font_Name,'FontSize',Label_Font_Size)
text(225,15,'Pr=0.71','Interpreter',Font_Interpreter,'Fontname',Font_Name,'FontSize',Label_Font_Size)
text(225,5,'Pr=0.10','Interpreter',Font_Interpreter,'Fontname',Font_Name,'FontSize',Label_Font_Size)

% plot the FDS results

H = 2;
dpdx = 9.0088e-6;
tau_w = 0.5*dpdx*H;
cp = 1;
mu = 1.8216e-5;
T_w = 20;

devcfile = {'heated_channel_Pr_0p10_16_devc.csv', ...
            'heated_channel_Pr_0p71_16_devc.csv', ...
            'heated_channel_Pr_1p00_16_devc.csv', ...
            'heated_channel_Pr_2p00_16_devc.csv'};

linefile = {'heated_channel_Pr_0p10_16_line.csv', ...
            'heated_channel_Pr_0p71_16_line.csv', ...
            'heated_channel_Pr_1p00_16_line.csv', ...
            'heated_channel_Pr_2p00_16_line.csv'};

skip_case = 0;

for i = [1,2,3,4]

    if ~exist([outdir,devcfile{i}])
        display(['Error: File ' [outdir,devcfile{i}] ' does not exist. Skipping case.'])
        skip_case = 1;
        continue
    end

    M = importdata([outdir,devcfile{i}],',',2);

    rho = M.data(end,strcmp(M.colheaders,'RHO')); % should be about 1.19
    u_tau = sqrt(tau_w/rho);
    delta_nu = (mu/rho)/u_tau;

    j1 = find(strcmp(M.colheaders,'HF0B'));
    j2 = find(strcmp(M.colheaders,'HF0T'));
    q_w = mean(M.data(end,j1:j2));
    T_tau = q_w/(rho*u_tau*cp);

    if ~exist([outdir,linefile{i}])
        display(['Error: File ' [outdir,linefile{i}] ' does not exist. Skipping case.'])
        skip_case = 1;
        continue
    end

    M = importdata([outdir,linefile{i}],',',2);

    zp = M.data(1:16,1)/delta_nu;

    switch(i)
        case {1,2,4}
            j1 = find(strcmp(M.colheaders,'T11'));
            j2 = find(strcmp(M.colheaders,'T115'));
            T1 = mean(M.data(1:16,j1:j2),2);      % bottom wall
            T2 = mean(M.data(32:-1:17,j1:j2),2);  % top wall
            Tp = (0.5*(T1+T2)-T_w)/T_tau;
            set(groot,'CurrentFigure',f2);
            hfig2(1)=semilogx(zp,Tp,'ksq-');
        case 3
            j1 = find(strcmp(M.colheaders,'U11'));
            j2 = find(strcmp(M.colheaders,'U115'));
            u1 = mean(M.data(1:16,j1:j2),2);     % bottom wall
            u2 = mean(M.data(32:-1:17,j1:j2),2); % top wall
            up = 0.5*(u1+u2)/u_tau;
            set(groot,'CurrentFigure',f1);
            hfig1(1)=semilogx(zp,up,'ksq-');
    end

    switch i
        case 1; err(i) = abs( mean(Tp)-mean(Tp_mean_Pr0p10,'omitnan') )/mean(Tp_mean_Pr0p10,'omitnan');
        case 2; err(i) = abs( mean(Tp)-mean(Tp_mean_Pr0p71,'omitnan') )/mean(Tp_mean_Pr0p71,'omitnan');
        case 3; err(i) = abs( mean(up)-mean(up_mean,'omitnan')        )/mean(up_mean,'omitnan');
        case 4; err(i) = abs( mean(Tp)-mean(Tp_mean_Pr2p00,'omitnan') )/mean(Tp_mean_Pr2p00,'omitnan');
    end

    if err(i)>1
        disp(['Matlab Warning: heated_channel case ',num2str(i),' out of tolerance'])
    end

end

if skip_case
    return
end

set(groot,'CurrentFigure',f1);
h = legend(hfig1,'FDS','DNS Re_\tau=180','Location','Northwest');
set(h,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [outdir,'heated_channel_Pr_1p00_16_git.txt'];
addverstr(a1,Git_Filename,'semilogx')

% print pdf

set(f1,'Visible',Figure_Visibility);
set(f1,'Units',Paper_Units)
set(f1,'PaperUnits',Paper_Units);
set(f1,'PaperSize',[Paper_Width Paper_Height]);
set(f1,'Position',[0 0 Paper_Width Paper_Height]);
print(f1,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/heated_channel_uplus')

set(groot,'CurrentFigure',f2);
h = legend(hfig2,'FDS','DNS Re_\tau=180','Location','Northwest');
set(h,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [outdir,'heated_channel_Pr_0p71_16_git.txt'];
addverstr(a2,Git_Filename,'semilogx')

% print pdf

set(f2,'Visible',Figure_Visibility);
set(f2,'Units',Paper_Units)
set(f2,'PaperUnits',Paper_Units);
set(f2,'PaperSize',[Paper_Width Paper_Height]);
set(f2,'Position',[0 0 Paper_Width Paper_Height]);
print(f2,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/heated_channel_Tplus')


% % Compute pressure gradient from Re_tau

% mu=1.8216e-5
% rho=1.1934
% Re_tau=180 % Kim and Moin, 1987
% delta=1 % channel half height
% u_tau=Re_tau*(mu/rho)/delta
% tau_w=rho*u_tau^2
% dpdx=tau_w/delta
% Pr = 1;
% Tref = 293;

% q = 2/(Pr*Re_tau) * rho*cp*u_tau*Tref/delta

% % Check bulk velocity

% figure
% for i = [1,2,3,4]
%     M = importdata([outdir,devcfile{i}],',',2);
%     t = M.data(:,find(strcmp(M.colheaders,'Time')));
%     UVEL = M.data(:,find(strcmp(M.colheaders,'"U-VEL"')));
%     plot(t,UVEL); hold on
% end

% % Check heat balance

% figure
% for i = [1,2,3,4]
%     M = importdata([outdir,devcfile{i}],',',2);
%     t = M.data(:,find(strcmp(M.colheaders,'Time')));
%     j1 = find(strcmp(M.colheaders,'"NetHF0B"'));
%     j2 = find(strcmp(M.colheaders,'"NetHF0T"'));
%     q_w = mean(M.data(:,j1:j2),2);
%     plot(t,q_w); hold on
% end





