%!/usr/bin/matlab
%Vanella from McDermott shunn_mms_error.m:
%03-09-2018

close all
clear all

% Shunn et al. Problem 3 parameters:
r0 = 5.;
r1 = 1.;
uf = 0.5;
vf = 0.5;
k  = 2.;
w  = 2.;
mu = 0.001;
D  = 0.001;

% Analytical solutions:
vd2d_mms_z = @(x,y,t) ...
( 1. + sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t) )/ ...
( (1+r0/r1) + (1-r0/r1)*sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t));

vd2d_mms_rho = @(x,y,t) ...
    1./( vd2d_mms_z(x,y,t)/r1 + (1-vd2d_mms_z(x,y,t))/r0 );

vd2d_mms_u = @(x,y,t) ...
    uf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*...
    cos(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*sin(pi*w*t);

vd2d_mms_v = @(x,y,t) ...
    vf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*...
    sin(pi*k*(x-uf*t))*cos(pi*k*(y-vf*t))*sin(pi*w*t);

vd2d_mms_H = @(x,y,t) ...
    0.5*(vd2d_mms_u(x,y,t)-uf)*(vd2d_mms_v(x,y,t)-vf);

vd2d_mms_p = @(x,y,t) ...
    vd2d_mms_rho(x,y,t)*(vd2d_mms_H(x,y,t) - ...
    0.5*(vd2d_mms_u(x,y,t)^2 + vd2d_mms_v(x,y,t)^2));

L = 2;
nx = [32,64,128,256,384];
dx = L./nx;

% % visualize field in time
% n=1;
% x=linspace(-L/2,L/2,nx(n)+1);
% xc = x(1:nx(n)) + .5*dx(n);
% yc = xc;
% for t=linspace(0,1,100)
%     for j=1:nx(n)
%         for i=1:nx(n)
%             %pres(i,j) = vd2d_mms_p(xc(i),yc(j),t);
%             hfld(i,j) = vd2d_mms_H(xc(i),yc(j),t);
%         end
%     end
%     surf(xc,yc,hfld)
%     axis([-1 1 -1 1 0 1])
%     pause(0.001)
% end
% return

datadir = '../../Verification/Complex_Geometry/';

file_s=1;
file_f=4;

file_str = ['exp_cen';'exp_gdv';'imp_cen';'imp_gdv'];

skip_case = 0;
for ifile=file_s:file_f
    
   file(ifile).name = { ['shunn3_32_cc_' file_str(ifile,:) '_mms.csv'], ...
                        ['shunn3_64_cc_' file_str(ifile,:) '_mms.csv'], ...
                        ['shunn3_128_cc_' file_str(ifile,:) '_mms.csv'],...
                        ['shunn3_256_cc_' file_str(ifile,:) '_mms.csv'],...
                        ['shunn3_384_cc_' file_str(ifile,:) '_mms.csv']};
   for n=1:length(file(ifile).name)
       if ~exist([datadir,file(ifile).name{n}])
           display(['Error: File ' [datadir,file(ifile).name{n}] ...
                    ' does not exist. Skipping case.'])
           skip_case = 1;
       end
   end
end
if skip_case
    return
end

for ifile=file_s:file_f
    
    e_r = zeros(length(file(ifile).name),1);
    e_z = zeros(length(file(ifile).name),1);
    e_u = zeros(length(file(ifile).name),1);
    e_H = zeros(length(file(ifile).name),1);
    
    
    for n=1:length(file(ifile).name)
        
        %disp(filename(n))
        
        M = importdata([datadir,file(ifile).name{n}],',',2); % FDS results.
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
        
        N = nx(n)*nx(n);
        
        e_r(n) = norm(e_r_vec(1:N),2)/nx(n);
        e_z(n) = norm(e_z_vec(1:N),2)/nx(n);
        e_u(n) = norm(e_u_vec(1:N),2)/nx(n);
        e_H(n) = norm(e_H_vec(1:N),2)/nx(n);
    end
    
    plot_style
    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    
    hh(1)=loglog(dx,e_r,'ko-');
    hold on
    hh(2)=loglog(dx,e_z,'ks-');
    hh(3)=loglog(dx,e_u,'k>-');
    hh(4)=loglog(dx,e_H,'k+-');
    hh(5)=loglog(dx,dx,'k--');
    hh(6)=loglog(dx,dx.^2,'k-');
    if (ifile==1 | ifile==3)
       axis([10^-3 10^-1 10^-6 10^-1])
    else
       axis([10^-3 10^-1 10^-6 10^0])
    end
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)
    
    xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',...
           Font_Interpreter,'Fontname',Font_Name)
    ylabel('{\it L_2} Error','FontSize',Title_Font_Size,'Interpreter', ...
           Font_Interpreter,'Fontname',Font_Name)
    lh=legend(hh,'FDS {\it \rho}','FDS {\it Z}','FDS {\it u}',...
              'FDS {\it H}','{\it O(\Deltax)}','{\it O(\Deltax^2)}',...
              'location','northwest');
    set(lh,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,...
        'Fontname',Font_Name)
    
    % add Git version if file is available
    Git_Filename = [datadir,'shunn3_256_cc_exp_cen_git.txt'];
    addverstr(gca,Git_Filename,'loglog')
    
    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    
    strng=['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/' ...
           'shunn_cc_' file_str(ifile,:) '_mms_convergence'];
    print(gcf,'-dpdf',strng)
    
    if (ifile==1 | ifile==3) % Central flux limiter cases:
        % check errors
        if e_r(end) > 2e-4
            display(['Matlab Warning: Density in shunn3 '    ...
            file_str(ifile,:) ' is out of tolerance. e_r = ',...
            num2str(e_r(end))])
        end
        if e_z(end) > 1e-4
            display(['Matlab Warning: Mixture fraction in shunn3 ' ...
            file_str(ifile,:) ' is out of tolerance. e_z = ',      ...
            num2str(e_z(end))])
        end
        if e_u(end) > 3e-5
            display(['Matlab Warning: Velocity in shunn3 '             ...
                     file_str(ifile,:) ' is out of tolerance. e_u = ', ...
                     num2str(e_u(end))])
        end
        if e_H(end) > 2e-3
            display(['Matlab Warning: Pressure in shunn3 '             ...
                     file_str(ifile,:) ' is out of tolerance. e_H = ', ...
                     num2str(e_H(end))])
        end
    end
    %disp([file_str(ifile,:) ', Err u=' num2str(e_u(end))])
end
