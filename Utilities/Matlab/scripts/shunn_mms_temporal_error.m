%!/usr/bin/matlab
%McDermott
%8-3-2015
%shunn_mms_temporal_error.m

close all
clear all

% Problem 1 parameters
r = 0.5;
nx = 256;

% Output files

datadir = '../../Verification/Scalar_Analytical_Solution/';
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
filename = {'shunn3_256_cfl_1_mms.csv', ...
            'shunn3_256_cfl_p5_mms.csv', ...
            'shunn3_256_cfl_p25_mms.csv', ...
            'shunn3_256_cfl_p125_mms.csv', ...
            'shunn3_256_cfl_p0625_mms.csv'};

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

% Gather FDS results

M1 = importdata([datadir,filename{3}],',',2);
M2 = importdata([datadir,filename{4}],',',2);
M3 = importdata([datadir,filename{5}],',',2);

rho_1 = M1.data(:,1);
rho_2 = M2.data(:,1);
rho_3 = M3.data(:,1);
Z_1 = M1.data(:,2);
Z_2 = M2.data(:,2);
Z_3 = M3.data(:,2);
U_1 = M1.data(:,3);
U_2 = M2.data(:,3);
U_3 = M3.data(:,3);
H_1 = M1.data(:,5);
H_2 = M2.data(:,5);
H_3 = M3.data(:,5);

p_rho = log( abs(rho_3-rho_2)./abs(rho_2-rho_1) )./log(r);
p_Z = log( abs(Z_3-Z_2)./abs(Z_2-Z_1) )./log(r);
p_U = log( abs(U_3-U_2)./abs(U_2-U_1) )./log(r);
p_H = log( abs(H_3-H_2)./abs(H_2-H_1) )./log(r);

% disp('Shunn 3 temporal order')
% disp(' ')
% disp(['L1 p rho = ',num2str( norm(p_rho,1)/(nx*nx) )])
% disp(['L2 p rho = ',num2str( norm(p_rho,2)/nx )])
% disp(['Linf p rho = ',num2str( norm(p_rho,-inf) )])
% disp(' ')
% disp(['L1 p Z = ',num2str( norm(p_Z,1)/(nx*nx) )])
% disp(['L2 p Z = ',num2str( norm(p_Z,2)/nx )])
% disp(['Linf p Z = ',num2str( norm(p_Z,-inf) )])
% disp(' ')
% disp(['L1 p U = ',num2str( norm(p_U,1)/(nx*nx) )])
% disp(['L2 p U = ',num2str( norm(p_U,2)/nx )])
% disp(['Linf p U = ',num2str( norm(p_U,-inf) )])
% disp(' ')
% disp(['L1 p H = ',num2str( norm(p_H,1)/(nx*nx) )])
% disp(['L2 p H = ',num2str( norm(p_H,2)/nx )])
% disp(['Linf p H = ',num2str( norm(p_H,-inf) )])
% disp(' ')

% flag errors

L2_rho = norm(p_rho,2)/nx;
if L2_rho<1.99
    disp(['Matlab Warning: L2_rho = ',num2str(L2_rho),' in Shunn 3 MMS temporal order'])
end

L2_Z = norm(p_Z,2)/sqrt(nx);
if L2_Z<1.99
    disp(['Matlab Warning: L2_Z = ',num2str(L2_Z),' in Shunn 3 MMS temporal order'])
end

% Tony Saad's way...
p1_rho_saad = log( norm(rho_3-rho_2,1)./norm(rho_2-rho_1,1) )./log(r);
p2_rho_saad = log( norm(rho_3-rho_2,2)./norm(rho_2-rho_1,2) )./log(r);
pinf_rho_saad = log( norm(rho_3-rho_2,inf)./norm(rho_2-rho_1,inf) )./log(r);
% disp(['L1 p rho Saad = ',num2str( p1_rho_saad )])
% disp(['L2 p rho Saad = ',num2str( p2_rho_saad )])
% disp(['Linf p rho Saad = ',num2str( pinf_rho_saad )])
% disp(' ')
p1_Z_saad = log( norm(Z_3-Z_2,1)./norm(Z_2-Z_1,1) )./log(r);
p2_Z_saad = log( norm(Z_3-Z_2,2)./norm(Z_2-Z_1,2) )./log(r);
pinf_Z_saad = log( norm(Z_3-Z_2,inf)./norm(Z_2-Z_1,inf) )./log(r);
% disp(['L1 p Z Saad = ',num2str( p1_Z_saad )])
% disp(['L2 p Z Saad = ',num2str( p2_Z_saad )])
% disp(['Linf p Z Saad = ',num2str( pinf_Z_saad )])
% disp(' ')
p1_U_saad = log( norm(U_3-U_2,1)./norm(U_2-U_1,1) )./log(r);
p2_U_saad = log( norm(U_3-U_2,2)./norm(U_2-U_1,2) )./log(r);
pinf_U_saad = log( norm(U_3-U_2,inf)./norm(U_2-U_1,inf) )./log(r);
% disp(['L1 p U Saad = ',num2str( p1_U_saad )])
% disp(['L2 p U Saad = ',num2str( p2_U_saad )])
% disp(['Linf p U Saad = ',num2str( pinf_U_saad )])
% disp(' ')
p1_H_saad = log( norm(H_3-H_2,1)./norm(H_2-H_1,1) )./log(r);
p2_H_saad = log( norm(H_3-H_2,2)./norm(H_2-H_1,2) )./log(r);
pinf_H_saad = log( norm(H_3-H_2,inf)./norm(H_2-H_1,inf) )./log(r);
% disp(['L1 p H Saad = ',num2str( p1_H_saad )])
% disp(['L2 p H Saad = ',num2str( p2_H_saad )])
% disp(['Linf p H Saad = ',num2str( pinf_H_saad )])

%surf(reshape(p_U,[nx,nx]),'EdgeColor','none','LineStyle','none','FaceLighting','phong')

% write the norms to a latex table

fid = fopen([plotdir,'shunn_terr_norms.tex'],'wt');

tmp_str = ['$L_1$ & ', ...
           num2str( norm(p_rho,1)/(nx*nx) ),' & ', ...
           num2str( norm(p_Z,1)/(nx*nx) ),' & ', ...
           num2str( norm(p_U,1)/(nx*nx) ),' & ', ...
           num2str( norm(p_H,1)/(nx*nx) ),' \\'];

fprintf(fid,'%s\n',tmp_str);

tmp_str = ['$L_2$ & ', ...
           num2str( norm(p_rho,2)/nx ),' & ', ...
           num2str( norm(p_Z,2)/nx ),' & ', ...
           num2str( norm(p_U,2)/nx ),' & ', ...
           num2str( norm(p_H,2)/nx ),' \\'];

fprintf(fid,'%s\n',tmp_str);

tmp_str = ['$L_1$ Saad & ', ...
           num2str( p1_rho_saad ),' & ', ...
           num2str( p1_Z_saad ),' & ', ...
           num2str( p1_U_saad ),' & ', ...
           num2str( p1_H_saad ),' \\'];

fprintf(fid,'%s\n',tmp_str);

tmp_str = ['$L_2$ Saad & ', ...
           num2str( p2_rho_saad ),' & ', ...
           num2str( p2_Z_saad ),' & ', ...
           num2str( p2_U_saad ),' & ', ...
           num2str( p2_H_saad ),' \\'];

fprintf(fid,'%s\n',tmp_str);

tmp_str = ['$L_{\infty}$ Saad & ', ...
           num2str( pinf_rho_saad ),' & ', ...
           num2str( pinf_Z_saad ),' & ', ...
           num2str( pinf_U_saad ),' & ', ...
           num2str( pinf_H_saad )];

fprintf(fid,'%s\n',tmp_str);

fclose(fid);







