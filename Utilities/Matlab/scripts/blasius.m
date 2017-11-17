% PARK HYUN WOOK, Yonsei University
% 8-15-2012
% blasius.m

close all
clear all

plot_style

repository = '../../Verification/Flowfields/';


%gather FDS result(n=16)
filename=[repository,'blasius_16_line.csv'];

if ~exist(filename)
    display(['Error: File ' filename ' does not exist. Skipping case.'])
    return
end

M = csvread(filename,2,0);

u_16 = M(:,2);
z_16 = M(:,1);


%gather FDS result(n=32)
filename=[repository,'blasius_32_line.csv'];

M = csvread(filename,2,0);

u_32 = M(:,2);
z_32 = M(:,1);


%gather FDS result(n=64)
filename=[repository,'blasius_64_line.csv'];

M = csvread(filename,2,0);

u_64 = M(:,2);
z_64 = M(:,1);

%gather blasius profile

u0= max(u_64(:));
zmax=0.3;
[eta,fp] = blasius_analytic(u0, zmax);

mu = 0.001;
rho = 1.19987607036;
xc = 0.05;
z_blasius=eta*sqrt(mu/rho*xc/u0);
u_blasius=fp*u0;


%plot whole velocity profile
range = 1:4:length(u_blasius);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(u_blasius(range),z_blasius(range),'bo');
hold on
H(2)=plot(u_16,z_16,'r--');
H(3)=plot(u_32,z_32,'c-.');
H(4)=plot(u_64,z_64,'g-');

axis([0 1.1 0 0.15])

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('{\it u} (m/s)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

h = legend(H,'Blasius','{\it N_z}=16','{\it N_z}=32','{\it N_z}=64','Location','northwest');
set(h,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

% add Git revision if file is available

Git_Filename = [repository,'blasius_16_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf for whole velocity profile
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_profile')


%%%%%%%get error comparing with analytic solution(blasius)%%%%%%%%%%
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
clear H

err(1) = 0;
err(2) = 0;
err(3) = 0;

% get error(n=16)
for i=1:length(u_16)
    err(1)=err(1)+(abs(u_16(i)-u_blasius(9+(i-1)*16)))^2;
end
err(1)=err(1)/16;
err(1)=sqrt(err(1));

% get error(n=32)
for i=1:length(u_32)
    err(2)=err(2)+(abs(u_32(i)-u_blasius(5+(i-1)*8)))^2 ;
end
err(2)=err(2)/32;
err(2)=sqrt(err(2));

% get error(n=64)
for i=1:length(u_64)
    err(3)=err(3)+(abs(u_64(i)-u_blasius(3+(i-1)*4)))^2;
end
err(3)=err(3)/64;
err(3)=sqrt(err(3));


dz(1)=abs(z_16(10)-z_16(9));
dz(2)=abs(z_32(10)-z_32(9));
dz(3)=abs(z_64(10)-z_64(9));

plot_style

H(1)=loglog(dz, err,'b*-','LineWidth',Line_Width); hold on
H(2)=loglog(dz, 10*dz,'k--','LineWidth',Line_Width);
H(3)=loglog(dz, 100*dz.^2,'k-','LineWidth',Line_Width);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Grid Spacing, {\it\deltaz} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('RMS Error (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([1e-3 1e-1 1e-3 1e-0])
legend_handle=legend(H,'FDS','{\itO}({\it\deltaz})','{\itO}({\it\deltaz^2})','Location','Northwest');
set(legend_handle,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available

Git_Filename = [repository,'blasius_16_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_convergence')



