% McDermott
% 12-7-2018
% soborot_mass_transport.m
%
% Solid body rotation flow field, scalar transport

close all
clear all

plot_style

data_dir = '../../Verification/Scalar_Analytical_Solution/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

chid = {'soborot_charm_cos_wave_16', ...
        'soborot_charm_cos_wave_32', ...
        'soborot_charm_cos_wave_64', ...
        'soborot_charm_square_wave_128', ...
        'soborot_charm_square_wave_16', ...
        'soborot_charm_square_wave_32', ...
        'soborot_charm_square_wave_64', ...
        'soborot_godunov_square_wave_128', ...
        'soborot_godunov_square_wave_16', ...
        'soborot_godunov_square_wave_32', ...
        'soborot_godunov_square_wave_64', ...
        'soborot_mp5_cos_wave_128', ...
        'soborot_mp5_cos_wave_16', ...
        'soborot_mp5_cos_wave_32', ...
        'soborot_mp5_cos_wave_64', ...
        'soborot_superbee_cos_wave_128', ...
        'soborot_superbee_cos_wave_16', ...
        'soborot_superbee_cos_wave_32', ...
        'soborot_superbee_cos_wave_64', ...
        'soborot_superbee_square_wave_128', ...
        'soborot_superbee_square_wave_128_1mesh', ...
        'soborot_superbee_square_wave_16', ...
        'soborot_superbee_square_wave_32', ...
        'soborot_superbee_square_wave_64'};


skip_case = 0;
for i=1:length(chid)
    devc_file = [data_dir,chid{i},'_devc.csv'];
    if ~exist(devc_file)
        display(['Error: File ' devc_file ' does not exist. Skipping case.'])
        skip_case = 1;
    end
end

if skip_case
    return
end

% positions of DEVC along diagonal D

L=1;
D=L*sqrt(2);
r_exact = linspace(0,D,1000);
Y_exact = zeros(1,length(r_exact));
i_range = find(r_exact>0.25 & r_exact<0.75);
Y_exact(i_range) = ones(1,length(i_range));
r_16 = (D/32):(D/16):(D-(D/32));
r_32 = (D/64):(D/32):(D-(D/64));
r_64 = (D/128):(D/64):(D-(D/128));
r_128 = (D/256):(D/128):(D-(D/256));

%%%%%% SQUARE WAVE %%%%%%

% FLUX_LIMITER='CHARM'

M = importdata([data_dir,'soborot_charm_square_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_charm_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_square_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_charm_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_square_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_charm_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_square_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_charm_128 = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_charm_16,'bo-');
H(3)=plot(r_32,Y_charm_32,'m*-');
H(4)=plot(r_64,Y_charm_64,'r^-');
H(5)=plot(r_128,Y_charm_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'CHARM','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_charm_square_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_charm_square_wave']);

% FLUX_LIMITER='SUPERBEE'

M = importdata([data_dir,'soborot_superbee_square_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_superbee_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_square_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_superbee_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_square_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_superbee_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_square_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_superbee_128 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_square_wave_128_1mesh_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_superbee_128_1mesh = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_superbee_16,'bo-');
H(3)=plot(r_32,Y_superbee_32,'m*-');
H(4)=plot(r_64,Y_superbee_64,'r^-');
H(5)=plot(r_128,Y_superbee_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'Superbee','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_superbee_square_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_superbee_square_wave']);

% FLUX_LIMITER='GODUNOV'

M = importdata([data_dir,'soborot_godunov_square_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_godunov_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_godunov_square_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_godunov_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_godunov_square_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_godunov_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_godunov_square_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_godunov_128 = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_godunov_16,'bo-');
H(3)=plot(r_32,Y_godunov_32,'m*-');
H(4)=plot(r_64,Y_godunov_64,'r^-');
H(5)=plot(r_128,Y_godunov_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'Godunov','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_godunov_square_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_godunov_square_wave']);

% plot error

Y_exact_16 = interp1(r_exact,Y_exact,r_16,'nearest');
Y_exact_32 = interp1(r_exact,Y_exact,r_32,'nearest');
Y_exact_64 = interp1(r_exact,Y_exact,r_64,'nearest');
Y_exact_128 = interp1(r_exact,Y_exact,r_128,'nearest');

% figure; plot(r_16,Y_exact_16,'ko'); hold on; plot([0.25 0.25],[0 1],'k-'); plot([0.75 0.75],[0 1],'k-')
% figure; plot(r_32,Y_exact_32,'ko'); hold on; plot([0.25 0.25],[0 1],'k-'); plot([0.75 0.75],[0 1],'k-')
% figure; plot(r_64,Y_exact_64,'ko'); hold on; plot([0.25 0.25],[0 1],'k-'); plot([0.75 0.75],[0 1],'k-')
% figure; plot(r_128,Y_exact_128,'ko'); hold on; plot([0.25 0.25],[0 1],'k-'); plot([0.75 0.75],[0 1],'k-')

e_charm_16 = norm(Y_charm_16-Y_exact_16,1)/length(r_16);
e_charm_32 = norm(Y_charm_32-Y_exact_32,1)/length(r_32);
e_charm_64 = norm(Y_charm_64-Y_exact_64,1)/length(r_64);
e_charm_128 = norm(Y_charm_128-Y_exact_128,1)/length(r_128);

e_superbee_16 = norm(Y_superbee_16-Y_exact_16,1)/length(r_16);
e_superbee_32 = norm(Y_superbee_32-Y_exact_32,1)/length(r_32);
e_superbee_64 = norm(Y_superbee_64-Y_exact_64,1)/length(r_64);
e_superbee_128 = norm(Y_superbee_128-Y_exact_128,1)/length(r_128);
e_superbee_128_1mesh = norm(Y_superbee_128_1mesh-Y_exact_128,1)/length(r_128);

% format long e
% abs(e_superbee_128_1mesh - e_superbee_128)
if abs(e_superbee_128_1mesh - e_superbee_128)>1e-10
    display(['Error: soborot_superbee_square_wave_128 single mesh and multi-mesh out of tolerance'])
end

e_godunov_16 = norm(Y_godunov_16-Y_exact_16,1)/length(r_16);
e_godunov_32 = norm(Y_godunov_32-Y_exact_32,1)/length(r_32);
e_godunov_64 = norm(Y_godunov_64-Y_exact_64,1)/length(r_64);
e_godunov_128 = norm(Y_godunov_128-Y_exact_128,1)/length(r_128);

dx = L./[16 32 64 128];

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=loglog(dx,dx.^.5,'k-.'); hold on
H(2)=loglog(dx,dx.^1,'k--');
H(3)=loglog(dx,[e_godunov_16 e_godunov_32 e_godunov_64 e_godunov_128],'ksq-');
H(4)=loglog(dx,[e_superbee_16 e_superbee_32 e_superbee_64 e_superbee_128],'k*-');
H(5)=loglog(dx,[e_charm_16 e_charm_32 e_charm_64 e_charm_128],'ko-');

axis([min(dx) .1 min(dx) max(dx.^.5)])
xlabel('Grid Spacing (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('L2 Error','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(.009,.18,'Square Wave Error','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'{\itO}(\delta{\itx}^{1/2})','{\itO}(\delta{\itx})','Godunov','Superbee','CHARM');
set(lh,'Location','Southeast')
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_superbee_square_wave_128_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_square_wave_error']);

%%%%%% COSINE WAVE %%%%%%

% FLUX_LIMITER='CHARM'

M = importdata([data_dir,'soborot_charm_cos_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_charm_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_cos_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_charm_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_cos_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_charm_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_charm_cos_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_charm_128 = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
Y_exact = zeros(1,length(r_exact));
i_range = find(r_exact>0.25 & r_exact<0.75);
Y_exact(i_range) = 0.5*(1. + cos(4.*pi*r_exact(i_range))); % see PERIODIC_TEST==13 in wall.f90
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_charm_16,'bo-');
H(3)=plot(r_32,Y_charm_32,'m*-');
H(4)=plot(r_64,Y_charm_64,'r^-');
H(5)=plot(r_128,Y_charm_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'CHARM','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_charm_cos_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_charm_cos_wave']);

% FLUX_LIMITER='SUPERBEE'

M = importdata([data_dir,'soborot_superbee_cos_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_superbee_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_cos_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_superbee_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_cos_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_superbee_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_superbee_cos_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_superbee_128 = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
Y_exact = zeros(1,length(r_exact));
i_range = find(r_exact>0.25 & r_exact<0.75);
Y_exact(i_range) = 0.5*(1. + cos(4.*pi*r_exact(i_range))); % see PERIODIC_TEST==13 in wall.f90
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_superbee_16,'bo-');
H(3)=plot(r_32,Y_superbee_32,'m*-');
H(4)=plot(r_64,Y_superbee_64,'r^-');
H(5)=plot(r_128,Y_superbee_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'Superbee','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_superbee_cos_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_superbee_cos_wave']);

% FLUX_LIMITER='MP5'

M = importdata([data_dir,'soborot_mp5_cos_wave_16_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-16"'));
Y_mp5_16 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_mp5_cos_wave_32_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-32"'));
Y_mp5_32 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_mp5_cos_wave_64_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-01"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-64"'));
Y_mp5_64 = M.data(end,col_start:col_end);

M = importdata([data_dir,'soborot_mp5_cos_wave_128_devc.csv'],',',2);
col_start = find(strcmp(M.colheaders,'"Y_TRACER-001"'));
col_end   = find(strcmp(M.colheaders,'"Y_TRACER-128"'));
Y_mp5_128 = M.data(end,col_start:col_end);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
Y_exact = zeros(1,length(r_exact));
i_range = find(r_exact>0.25 & r_exact<0.75);
Y_exact(i_range) = 0.5*(1. + cos(4.*pi*r_exact(i_range))); % see PERIODIC_TEST==13 in wall.f90
H(1)=plot(r_exact,Y_exact,'k-'); hold on
H(2)=plot(r_16,Y_mp5_16,'bo-');
H(3)=plot(r_32,Y_mp5_32,'m*-');
H(4)=plot(r_64,Y_mp5_64,'r^-');
H(5)=plot(r_128,Y_mp5_128,'gsq-');

xlabel('Radial Position (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Scalar Mass Fraction','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([0 1 0 1.2])
text(.05,1.1,'MP5','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exact','{\itn}=16','{\itn}=32','{\itn}=64','{\itn}=128');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_mp5_cos_wave_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_mp5_cos_wave']);

% plot error

Y_exact_16 = interp1(r_exact,Y_exact,r_16);
Y_exact_32 = interp1(r_exact,Y_exact,r_32);
Y_exact_64 = interp1(r_exact,Y_exact,r_64);
Y_exact_128 = interp1(r_exact,Y_exact,r_128);

e_charm_16 = norm(Y_charm_16-Y_exact_16,1)/length(r_16);
e_charm_32 = norm(Y_charm_32-Y_exact_32,1)/length(r_32);
e_charm_64 = norm(Y_charm_64-Y_exact_64,1)/length(r_64);
e_charm_128 = norm(Y_charm_128-Y_exact_128,1)/length(r_128);

e_superbee_16 = norm(Y_superbee_16-Y_exact_16,1)/length(r_16);
e_superbee_32 = norm(Y_superbee_32-Y_exact_32,1)/length(r_32);
e_superbee_64 = norm(Y_superbee_64-Y_exact_64,1)/length(r_64);
e_superbee_128 = norm(Y_superbee_128-Y_exact_128,1)/length(r_128);

e_mp5_16 = norm(Y_mp5_16-Y_exact_16,1)/length(r_16);
e_mp5_32 = norm(Y_mp5_32-Y_exact_32,1)/length(r_32);
e_mp5_64 = norm(Y_mp5_64-Y_exact_64,1)/length(r_64);
e_mp5_128 = norm(Y_mp5_128-Y_exact_128,1)/length(r_128);

dx = L./[16 32 64 128];

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=loglog(dx,dx.^1,'k-.'); hold on
H(2)=loglog(dx,dx.^2,'k--');
H(3)=loglog(dx,[e_charm_16 e_charm_32 e_charm_64 e_charm_128],'ko-');
H(4)=loglog(dx,[e_superbee_16 e_superbee_32 e_superbee_64 e_superbee_128],'k*-');
H(5)=loglog(dx,[e_mp5_16 e_mp5_32 e_mp5_64 e_mp5_128],'ksq-');

axis([min(dx) .1 min(dx.^2) max(dx)])
xlabel('Grid Spacing (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('L2 Error','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(.009,.04,'Cosine Wave Error','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'{\itO}(\delta{\itx})','{\itO}(\delta{\itx}^2)','CHARM','Superbee','MP5');
set(lh,'Location','Southeast')
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
legend boxoff

Git_Filename = [data_dir,'soborot_superbee_cos_wave_128_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'soborot_cos_wave_error']);

