% McDermott
% 7-7-14
% mccaffrey_plumes.m

close all
clear all

plot_style
Font_Interpreter='latex';

Q = [14.4 21.7 33.0 44.9 57.5]; % kW [14.4 21.7 33.0 44.9 57.5]
g = 9.8;
rho = 1.18;
cp = 1;
T0 = 273.15 + 20;

DS = (Q/(rho*cp*T0*sqrt(g))).^(2/5); % m

datadir = '../../Validation/McCaffrey_Plume/FDS_Output_Files/';
plotdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/McCaffrey_Plume/';

%chid = {'McCaffrey_14_kW_10','McCaffrey_22_kW_10','McCaffrey_33_kW_10','McCaffrey_45_kW_10','McCaffrey_57_kW_10'};
%chid = {'McCaffrey_14_kW_20','McCaffrey_22_kW_20','McCaffrey_33_kW_20','McCaffrey_45_kW_20','McCaffrey_57_kW_20'};
chid = {'McCaffrey_14_kW_40','McCaffrey_22_kW_40','McCaffrey_33_kW_40','McCaffrey_45_kW_40','McCaffrey_57_kW_40'};
mark = {'ko','k+','k^','ksq','kd'};
n_chid = length(chid);

% McCaffrey plume correlations

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

% Baum and McCaffrey plume correlations (in terms of D*)
% 2nd IAFSS, pp. 129-148

zs = logspace(-2,1,100);
n = [1/2 0 -1/3];
A = [2.18 2.45 3.64];
B = [2.91 3.81 8.41];

for i=1:length(zs)
    if zs(i)<1.32
        j=1;
    elseif zs(i)>=1.32 & zs(i)<=3.3
        j=2;
    elseif zs(i)>3.3
        j=3;
    end
    us(i) = A(j)*zs(i)^n(j);
    Ts(i) = B(j)*zs(i)^(2*n(j)-1);
end

hh(n_chid+1)=loglog(zs,us,'b--','linewidth',2); hold on
xmin = 0.2;
xmax = 20;
ymin = 1;
ymax = 3.5;
axis([xmin xmax ymin ymax])
%grid on
%xlabel('$z/Q^{2/5}$ (m $\cdot$ kW$^{-2/5}$ )','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
%ylabel('$V/Q^{1/5}$ (m $\cdot$ s$^{-1}$ $\cdot$ kW$^{-1/5}$ )','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z^* = z/D^*$','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('$v^* = v/\sqrt{g D^*} $','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

% FDS results velocity

for i=1:n_chid
    M = importdata([datadir,chid{i},'_line.csv'],',',2);
    z = M.data(:,1);
    v = M.data(:,3);

    %zq_fds = z./Q(i)^(2/5);
    %vq_fds = v./Q(i)^(1/5);

    zs_fds = z./DS(i);
    vs_fds = v./sqrt(g*DS(i));

    hh(i)=loglog(zs_fds,vs_fds,mark{i});
end

leg_key = {'14.4 kW','21.7 kW','33.0 kW','44.9 kW','57.5 kW','$A(z^*\,)^n$'};
%leg_key = {'57.5 kW','$A(z^*\,)^n$'};
lh = legend(hh,leg_key,'location','northwest');
set(lh,'Interpreter',Font_Interpreter)
legend boxoff

% add SVN if file is available

SVN_Filename = [datadir,'McCaffrey_14_kW_10_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plotdir,'McCaffrey_Velocity_Correlation'])


% FDS results temperature rise

figure

hh(n_chid+1)=loglog(zs,Ts,'r--','linewidth',2); hold on

xmin = 0.1;
xmax = 20;
ymin = 0.1;
ymax = 5;
axis([xmin xmax ymin ymax])
%grid on
%xlabel('$z/Q^{2/5}$ (m $\cdot$ kW$^{-2/5}$ )','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
%ylabel('$\Delta T$ ($^\circ$C)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
xlabel('$z^* = z/D^*$','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
ylabel('$\Theta^* = (T-T_0)/T_0$','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

for i=1:n_chid
    M = importdata([datadir,chid{i},'_line.csv'],',',2);
    z = M.data(:,1);
    T = M.data(:,2) + 273.15;
    %zq_fds = z./Q(i)^(2/5);
    %hh(i)=loglog(zq_fds,T,mark{i});
    hh(i)=loglog(z./DS(i),(T-T0)/T0,mark{i});
end

leg_key = {'14.4 kW','21.7 kW','33.0 kW','44.9 kW','57.5 kW','$B(z^*\,)^{2n-1}$'};
%leg_key = {'57.5 kW','$B(z^*\,)^{2n-1}$'};
lh = legend(hh,leg_key,'location','southwest');
set(lh,'Interpreter',Font_Interpreter)
legend boxoff

% add SVN if file is available

SVN_Filename = [datadir,'McCaffrey_14_kW_10_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plotdir,'McCaffrey_Temperature_Correlation'])




