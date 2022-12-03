% McDermott and Henderson
% 26 Oct 2017 (updated 11-27-2022 RJM)
% natconh.m

close all
clear all

plot_style

results_dir = ['../../../out/Convection/'];

% calculations below were used for input file setup

g = 9.80665;
%       1      2     3     4     5     6     7    8    9   10   11   12   13   14   15   16   17   18   19
S = [.002   .002   .02   .02   .02   .02   .02   .2   .2   .2   .2    2    2    2    2   20   20   20   20];
T1 = [303    600   295   297   303   333   503  295  303  333  503  295  303  333  503  295  303  333  503];
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

Ra = (g*beta.*(T1-T2).*S.^3)./(alpha.*nu);

% see J.P. Holman p. 361 for correlations

Ra_Limit_1=1183;
Ra_Limit_2=5049;
Ra_Limit_3=3.10508e6;

r0=[];
r1=[];
r2=[];
r3=[];

for i=1:length(Ra)
    if Ra(i)<Ra_Limit_1
        r0=[r0 i];
        Nu(i)=1;
    elseif Ra(i)>Ra_Limit_1 & Ra(i)<=Ra_Limit_2
        r1=[r1 i];
        Nu(i)=0.059*Ra(i)^(.4);
    elseif Ra(i)>Ra_Limit_2 & Ra(i)<=Ra_Limit_3
        r2=[r2 i];
        Nu(i)=0.212*Ra(i)^(.25);
    elseif Ra(i)>Ra_Limit_3
        r3=[r3 i];
        Nu(i)=0.061*Ra(i)^(1/3);
    end
end

RAYLEIGH = logspace(0,14,1000);
for i=1:length(RAYLEIGH)
    if RAYLEIGH(i)<Ra_Limit_1
        NUSSELT(i)=1;
    elseif RAYLEIGH(i)>Ra_Limit_1 & RAYLEIGH(i)<=Ra_Limit_2
        NUSSELT(i)=0.059*RAYLEIGH(i)^(.4);
    elseif RAYLEIGH(i)>Ra_Limit_2 & RAYLEIGH(i)<=Ra_Limit_3
        NUSSELT(i)=0.212*RAYLEIGH(i)^(.25);
    elseif RAYLEIGH(i)>Ra_Limit_3
        NUSSELT(i)=0.061*RAYLEIGH(i)^(1/3);
    end
end

figure
set(gcf,'Visible',Figure_Visibility);
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% loglog(Ra(r0),Nu(r0),'k-'); hold on
% loglog(Ra(r1),Nu(r1),'k-')
% loglog(Ra(r2),Nu(r2),'k-')
% loglog(Ra(r3),Nu(r3),'k-')
marker_handle(1)=loglog(RAYLEIGH,NUSSELT,'k-'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
axis([1 1e15 .5 1e4])
xlabel('Rayleigh Number','FontSize',Label_Font_Size)
ylabel('Nusselt Number','FontSize',Label_Font_Size)
%title('Natural Convection in Horizontal Enclosures')

% FDS results

casename={...
'natconh_1',...
'natconh_2',...
'natconh_3',...
'natconh_4',...
'natconh_5',...
'natconh_6',...
'natconh_7',...
'natconh_8',...
'natconh_9',...
'natconh_10',...
'natconh_11',...
'natconh_12',...
'natconh_13',...
'natconh_14',...
'natconh_15',...
'natconh_16',...
'natconh_17',...
'natconh_18',...
'natconh_19'...
};

delta = S;
T     = T1;

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

        Q = mean(mean(abs(M.data(:,col))*1000));  % heat flow, W/m2

        Nu_FDS = Q*(delta(i)/k)/(T(i)-T2);

        marker_handle(j+1)=plot(Ra(i),Nu_FDS,marker_style{j},'MarkerSize',8);
    end
end

lh=legend(marker_handle,'Correlation, Nu = C Ra^n','FDS {\itH/\Deltax}=8','FDS {\itH/\Deltax}=16','FDS {\itH/\Deltax}=32','Location','Northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'natconh_19_32_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natural_convection_horizontal_enclosure');


% Complex geometry cases (rotated grid 18 degrees)

if 1

figure
set(gcf,'Visible',Figure_Visibility);
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% loglog(Ra(r0),Nu(r0),'k-'); hold on
% loglog(Ra(r1),Nu(r1),'k-')
% loglog(Ra(r2),Nu(r2),'k-')
% loglog(Ra(r3),Nu(r3),'k-')
marker_handle(1)=loglog(RAYLEIGH,NUSSELT,'k-'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
axis([1 1e15 .5 1e4])
xlabel('Rayleigh Number','FontSize',Label_Font_Size)
ylabel('Nusselt Number','FontSize',Label_Font_Size)
%title('Natural Convection in Horizontal Enclosures')

% FDS results

casename={...
'natconh_1',...
'natconh_2',...
'natconh_3',...
'natconh_4',...
'natconh_5',...
'natconh_6',...
'natconh_7',...
'natconh_8',...
'natconh_9',...
'natconh_10',...
'natconh_11',...
'natconh_12',...
'natconh_13',...
'natconh_14',...
'natconh_15',...
'natconh_16',...
'natconh_17',...
'natconh_18',...
'natconh_19'...
};

delta = S;
T     = T1;

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
        A = (delta(i)*8)^2;
        q = mean(abs(Q2(qrange)))*1000/A;  % heat flux, W/m2

        Nu_FDS = q*(delta(i)/k)/(T(i)-T2);

        marker_handle(j+1)=plot(Ra(i),Nu_FDS,marker_style{j},'MarkerSize',8);
    end
end

lh=legend(marker_handle(1:3),'Correlation, Nu = C Ra^n','GEOM {\itH/\Deltax}=8','GEOM {\itH/\Deltax}=16','Location','Northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'natconh_19_16_rot_18_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natconh_geom');

end
