% PARK HYUN WOOK, Yonsei University
% 8-15-2012
% blasius.m

close all
clear all

repository = '../../../Verification/Flowfields/';


%gather FDS result(n=16)
filename=[repository,'blasius_16_line.csv'];

M = csvread(filename,2,0);

u_16 = M(:,2);
z_16 = M(:,1);


%gather FDS result(n=32)
filename=[repository,'blasius_32_line.csv'];

M = csvread(filename,2,0);

u_32 = M(:,2);
z_32 = M(:,1);


%gather FDS result(n=64)
filename=[repository,'blasius_32_line.csv'];

M = csvread(filename,2,0);

u_64 = M(:,2);
z_64 = M(:,1);


%gather FDS result(n=128)
filename=[repository,'blasius_32_line.csv'];

M = csvread(filename,2,0);

u_128 = M(:,2);
z_128 = M(:,1);


% %gather FDS result(u0)
% filename=[repository,'blasius_128_devc.csv'];
% 
% M = csvread(filename,2,0);


%gather blasius profile

u0= max(u_128(:));
zmax=0.3;
[eta,fp] = blasius_analytic(u0, zmax);

z_blasius=eta*sqrt(1E-3/1.199*0.05/u0);
u_blasius=fp*u0;


%plot whole velocity profile
range = 1:4:length(u_blasius);
H(1)=plot(u_blasius(range),z_blasius(range),'bo');
hold on

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

axis([0 1.1 0 0.15]) 
xlabel('$u$','Interpreter',Font_Interpreter,'FontSize',Title_Font_Size)
ylabel('$z$','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'Rotation',0.0)


H(2)=plot(u_16,z_16,'-.r');
H(3)=plot(u_32,z_32,'--c');
H(4)=plot(u_64,z_64,'-.g');
H(5)=plot(u_128,z_128,'-k');


h = legend(H,'Blasius','$N_z=16$','$N_z=32$','$N_z=64$','$N_z=128$','Location','northwest');
legend boxoff
set(h,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% print to pdf for whole velocity profile
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[1.1*Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 1.1*Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_profile')


%%%%%%%get error comparing with analytic solution(blasius)%%%%%%%%%%
figure
clear H

err(1) = 0;
err(2) = 0;
err(3) = 0;
err(4) = 0;

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

% get error(n=128)
for i=1:length(u_128)   
    err(4)=err(4)+(abs(u_128(i)-u_blasius(2+(i-1)*2)))^2;    
end   
err(4)=err(4)/length(u_128);
err(4)=sqrt(err(4));


dz(1)=abs(z_16(10)-z_16(9));
dz(2)=abs(z_32(10)-z_32(9));
dz(3)=abs(z_64(10)-z_64(9));
dz(4)=abs(z_128(10)-z_128(9));

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

H(1)=loglog(dz, err,'b*-','LineWidth',Line_Width); hold on
H(2)=loglog(dz, dz,'k--','LineWidth',Line_Width);
H(3)=loglog(dz, dz.^2,'k-','LineWidth',Line_Width);

xlabel('Grid Spacing, $\delta z$ (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('L2 Error (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
%axis([1e-3 1e-1 1e-3 1e-0])
legend_handle=legend(H,'FDS','${\cal O}(\delta z)$','${\cal O}(\delta z^2)$','Location','Northwest');
legend boxoff
set(legend_handle,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[1.1*Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 1.1*Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_convergence')



