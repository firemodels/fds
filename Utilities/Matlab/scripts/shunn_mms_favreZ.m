%!/usr/bin/matlab
%McDermott
%09-04-2013

close all
clear all

% Shunn et al. Problem 3 parameters
r0 = 5.;
r1 = 1.;
uf = 0.5;
vf = 0.5;
k = 2.;
w = 2.;
mu = 0.001;
D  = 0.001;

% analytical solutions
vd2d_mms_z = @(x,y,t) ...
    ( 1. + sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t) )/ ...
    ( (1+r0/r1) + (1-r0/r1)*sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t) );

vd2d_mms_rho = @(x,y,t) ...
    1./( vd2d_mms_z(x,y,t)/r1 + (1-vd2d_mms_z(x,y,t))/r0 );

X = 0;
Y = 0;
T = 10;
N = 1000;
dt = T/N;

Zbar    = 0;
rhoZbar = 0;
rhobar  = 0;
for n=1:N
    t = n*dt;
    Zbar    = Zbar    + vd2d_mms_z(X,Y,t)*dt;
    rhoZbar = rhoZbar + vd2d_mms_rho(X,Y,t)*vd2d_mms_z(X,Y,t)*dt;
    rhobar  = rhobar  + vd2d_mms_rho(X,Y,t)*dt;
end

Zbar_mms = Zbar/T;
rhoZbar = rhoZbar/T;
rhobar = rhobar/T;

FavreZ_mms = rhoZbar/rhobar;

% read the FDS output

datadir = '../../Verification/Scalar_Analytical_Solution/';
filename = {'shunn3_FavreZ_32_devc.csv','shunn3_FavreZ_64_devc.csv'};
linestyle = {'k-','k--','m-','m--','b-','b--'};

skip_case = 0;

for n=1:length(filename)
    if ~exist([datadir,filename{n}])
        display(['Error: File ' [datadir,finlename{n}] ' does not exist. Skipping case.'])
        skip_case = 1;
    end
end

if skip_case
    return
end

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hh(1)=plot([0,T],[FavreZ_mms,FavreZ_mms],'k-','linewidth',2); hold on
hh(2)=plot([0,T],[Zbar_mms,Zbar_mms],'k--','linewidth',2);

k=2;
for n=1:length(filename)
    M = importdata([datadir,filename{n}],',',2);
    t = M.data(:,1);
    Zbar = M.data(:,2);
    FavreZ = M.data(:,3);
    k=k+1;
    hh(k) = plot(t,FavreZ,linestyle{k});
    k=k+1;
    hh(k) = plot(t,Zbar,linestyle{k});
end

axis([0 T 0 0.15])

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('Mass Fraction','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
lh=legend(hh,'Analytical Favre','Analytical Average', ...
    'FDS 32 Favre','FDS 32 Average','FDS 64 Favre','FDS 64 Average','location','southeast');
set(lh,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)

% add Git version if file is available

Git_Filename = [datadir,'shunn3_FavreZ_32_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/shunn_mms_FavreZ')

% check errors
FavreZ_error_64 = abs(FavreZ(end)-FavreZ_mms)/FavreZ_mms;
if FavreZ_error_64 > 0.01
   display(['Matlab Warning: FavreZ in shunn3_FavreZ_64 is out of tolerance. Error = ',num2str(FavreZ_error_64)])
end


