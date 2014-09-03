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

vd2d_mms_u = @(x,y,t) ...
    uf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*cos(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*sin(pi*w*t);

vd2d_mms_v = @(x,y,t) ...
    vf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*sin(pi*k*(x-uf*t))*cos(pi*k*(y-vf*t))*sin(pi*w*t);

vd2d_mms_p = @(x,y,t) ...
    .5*vd2d_mms_rho(x,y,t)*vd2d_mms_u(x,y,t)*vd2d_mms_v(x,y,t);

vd2d_mms_H = @(x,y,t) ...
    vd2d_mms_p(x,y,t)/vd2d_mms_rho(x,y,t) + 0.5*(vd2d_mms_u(x,y,t)^2+vd2d_mms_v(x,y,t)^2);

L = 2;
nx = [32,64,128,256,512];
dx = L./nx;

datadir = '../../Verification/Scalar_Analytical_Solution/';
filename = {'shunn3_32_mms.csv','shunn3_64_mms.csv','shunn3_128_mms.csv','shunn3_256_mms.csv','shunn3_512_mms.csv'};

skip_case = 0;

if ~exist([datadir,'shunn3_32_mms.csv'])
    display(['Error: File ' [datadir,'shunn3_32_mms.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([datadir,'shunn3_64_mms.csv'])
    display(['Error: File ' [datadir,'shunn3_64_mms.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([datadir,'shunn3_128_mms.csv'])
    display(['Error: File ' [datadir,'shunn3_128_mms.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([datadir,'shunn3_256_mms.csv'])
    display(['Error: File ' [datadir,'shunn3_256_mms.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([datadir,'shunn3_512_mms.csv'])
    display(['Error: File ' [datadir,'shunn3_512_mms.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

e_r = zeros(length(filename),1);
e_z = zeros(length(filename),1);
e_u = zeros(length(filename),1);
e_H = zeros(length(filename),1);


for n=1:length(filename)

    %disp(filename(n))

    M = importdata([datadir,filename{n}],',',2); % collect FDS results
    T = str2num(M.textdata{2,1}); % exact time value

    x=linspace(-L/2,L/2,nx(n)+1);
    xc = x(1:nx(n)) + .5*dx(n);
    yc = xc;

    % inialize error arrays
    rho_error = zeros(nx(n),nx(n));
    e_r_vec = zeros(nx(n)*nx(n),1);

    z_error = zeros(nx(n),nx(n));
    e_z_vec = zeros(nx(n)*nx(n),1);

    u_error = zeros(nx(n),nx(n));
    e_u_vec = zeros(nx(n)*nx(n),1);

    H_error = zeros(nx(n),nx(n));
    e_H_vec = zeros(nx(n)*nx(n),1);

    p = 0;
    for j=1:nx(n)
        for i=1:nx(n)
            p = p+1;
            
            rho = M.data(p,1);
            rho_mms = vd2d_mms_rho(xc(i),yc(j),T);
            rho_error(i,j) = rho - rho_mms;
            e_r_vec(p) = rho_error(i,j);
            
            z = M.data(p,2);
            z_mms = vd2d_mms_z(xc(i),yc(j),T);
            z_error(i,j) = z - z_mms;
            e_z_vec(p) = z_error(i,j);
            
            u = M.data(p,3);
            u_mms = vd2d_mms_u(x(i+1),yc(j),T);
            u_error(i,j) = u - u_mms;
            e_u_vec(p) = u_error(i,j);
            
            H = M.data(p,5);
            H_mms = vd2d_mms_H(xc(i),yc(j),T);
            H_error(i,j) = H - H_mms;
            e_H_vec(p) = H_error(i,j);
        end
    end

    e_r(n) = norm(e_r_vec,2)/(nx(n)^2);
    e_z(n) = norm(e_z_vec,2)/(nx(n)^2);
    e_u(n) = norm(e_u_vec,2)/(nx(n)^2);
    e_H(n) = norm(e_H_vec,2)/(nx(n)^2);
end

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

hh(1)=loglog(dx,e_r,'ko-');
hold on
hh(2)=loglog(dx,e_z,'ks-');
hh(3)=loglog(dx,e_u,'k>-');
hh(4)=loglog(dx,e_H,'k+-');
hh(5)=loglog(dx,0.1*dx,'k--');
hh(6)=loglog(dx,dx.^2,'k-');
axis([10^-3 10^-1 10^-7 10^-1])

xlabel('{\it \Delta x} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('L2 Error','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
legend(hh,'FDS {\it \rho}','FDS {\it Z}','FDS {\it u}','FDS {\it H}','{\it O(\Delta x)}','{\it O(\Delta x^2)}','location','northwest')
legend('boxoff')

% add SVN if file is available

SVN_Filename = [datadir,'shunn3_256_svn.txt'];
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
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/shunn_mms_convergence')

% check errors
if e_r(length(e_r)) > 4.96e-06
   display(['Matlab Warning: Density in shunn3 is out of tolerance.'])
end
if e_z(length(e_z)) > 9.57e-07
   display(['Matlab Warning: Mixture fraction in shunn3 is out of tolerance.'])
end
if e_u(length(e_u)) > 3.34e-06
   display(['Matlab Warning: Velocity in shunn3 is out of tolerance.'])
end
if e_H(length(e_H)) > 7.63e-04
   display(['Matlab Warning: Pressure in shunn3 is out of tolerance.'])
end


