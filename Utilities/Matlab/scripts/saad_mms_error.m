%!/usr/bin/matlab
%McDermott
%7-14-2015
%saad_mms_error.m

close all
clear all

% Problem 1 parameters
r0 = 5;
r1 = .5;
u0 = 1;
sig = 0.1;

L = 2;
nx = [32,64,128,256,512,1024];
ny = [4,4,4,4,4,4];
dx = L./nx;

% analytical solutions

%f0 = @(x) exp(-.5*x^2/sig^2)
f0 = @(x) 0.5*(1+sin(2*pi*x/L))

vd1d_mms_z = @(x,t,u0) f0(x-u0*t)

vd1d_mms_rho = @(x,t,u0) 1/( vd1d_mms_z(x,t,u0)/r1 + (1-vd1d_mms_z(x,t,u0))/r0 )

datadir = '../../Verification/Scalar_Analytical_Solution/';
filename = {'saad_32_mms.csv','saad_64_mms.csv','saad_128_mms.csv','saad_256_mms.csv','saad_512_mms.csv','saad_1024_mms.csv'};

skip_case = 0;

for n=1:length(filename)
    if ~exist([datadir,filename{n}])
        display(['Error: File ' [datadir,filename{n}] ' does not exist. Skipping case.'])
        skip_case = 1;
    end
end

if skip_case
    return
end

e_r = zeros(length(filename),1);
e_z = zeros(length(filename),1);

for n=1:length(filename)

    %disp(filename(n))

    M = importdata([datadir,filename{n}],',',2); % collect FDS results
    T = str2num(M.textdata{2,1}); % exact time value
    U(n) = u0; %mean(M.data(:,3));

    x=linspace(-L/2,L/2,nx(n)+1);
    xc = x(1:nx(n)) + .5*dx(n);

    % inialize error arrays
    rho_error = zeros(nx(n),ny(n));
    e_r_vec = zeros(nx(n)*ny(n),1);

    z_error = zeros(nx(n),ny(n));
    e_z_vec = zeros(nx(n)*ny(n),1);

    p = 0;
    for j=1:ny(n)
        for i=1:nx(n)
            p = p+1;

            rho = M.data(p,1);
            rho_mms = vd1d_mms_rho(xc(i),T,U(n));
            rho_error(i,j) = rho - rho_mms;
            e_r_vec(p) = rho_error(i,j);

            z = M.data(p,2);
            z_mms = vd1d_mms_z(xc(i),T,U(n));
            z_error(i,j) = z - z_mms;
            e_z_vec(p) = z_error(i,j);
        end
    end

    N = nx(n)*ny(n);
    e_r(n) = norm(e_r_vec,2)/sqrt(N);
    e_z(n) = norm(e_z_vec,2)/sqrt(N);
end

sqrt(e_z(end-1)/e_z(end))
sqrt(e_r(end-1)/e_r(end))

plot_style
figure
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

hh(1)=loglog(dx,e_r,'ko-');
hold on
hh(2)=loglog(dx,e_z,'ks-');
hh(3)=loglog(dx,dx,'k--');
hh(4)=loglog(dx,dx.^2,'k-');
%axis([10^-3 10^-1 10^-7 10^-1])

xlabel('{\it \Delta x} (m)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('L2 Error','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
lh=legend(hh,'FDS {\it \rho}','FDS {\it Z}','{\it O(\Delta x)}','{\it O(\Delta x^2)}','location','northwest')
set(lh,'FontSize',Key_Font_Size)
set(lh,'Interpreter',Font_Interpreter)

% add version string if file is available

Git_Filename = [datadir,'saad_256_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/saad_mms_convergence')

% check errors
if e_r(length(e_r)) > 1
   display(['Matlab Warning: Density in saad is out of tolerance. e_r = ',num2str(e_r(length(e_r)))])
end
if e_z(length(e_z)) > 1
   display(['Matlab Warning: Mixture fraction in saad is out of tolerance. e_z = ',num2str(e_z(length(e_z)))])
end


