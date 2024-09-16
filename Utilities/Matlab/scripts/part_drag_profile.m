% McDermott
% 8-27-24
% part_drag_profile.m

close all
clear all

plot_style

figure
set(gcf,'Visible',Figure_Visibility);
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% prescribed velocity profile

% mpv = 4.;                 % kg/m^3, mass_per_volume (from FDS input file)
% v_xb = 10^3;              % volume of XB region on init line in FDS input file
% nppc = 10;                % number of particles per cell
% n = 5*5*20*nppc;          % number of particles
% rho_p = 400;              % density of grass, kg/m^3
r_p = 0.001;              % radius, m
l_p = 0.02;               % length, m
v_p = pi*(r_p)^2*l_p;     % volume of a single particle, m^3
shape_factor = 0.25;                % assumes random orientation of cylinders
a_p = shape_factor*l_p*(2*pi*r_p);  % projected area, m^2
% m_p = rho_p*v_p;          % mass of single particle, kg
% pwt = mpv*v_xb/(n*m_p)    % particle weight factor

z   = linspace(0,10,20);
u_z = z;
c_d = 2.8;       % from FDS input file (specified)
rho_g = 1.195;   % from FDS out file
f_x = c_d * a_p * 0.5*rho_g*(u_z.^2);  % drag experienced by a single particle

H(1)=plot(z,-f_x,'k-'); hold on

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

ddir='../../Verification/WUI/';
chid={'part_drag_prof_ux','part_drag_prof_uy','part_drag_prof_uz',...
      'part_drag_prof_vx','part_drag_prof_vy','part_drag_prof_vz',...
      'part_drag_prof_wx','part_drag_prof_wy','part_drag_prof_wz'};
j={1,2,3,1,2,3,1,2,3}; % coordinate direction (x=1, y=2, z=3)

for i=1:length(chid) % chid_for

    skip_case = 0;
    if ~exist([ddir,chid{i},'_1.prt5'])
        display(['Error: File ' [ddir,chid{i},'_1.prt5'] ' does not exist. Skipping case.'])
        skip_case = 1;
    end

    if skip_case
        return
    end

    [STIME, XP, YP, ZP, QP] = read_prt5([ddir,chid{i},'_1.prt5'],'real*4');

    switch j{i}
        case 1
            H(2)=plot(XP(end,:),QP(end,:,1,1)./QP(end,:,1,2),'b.');
            v = abs( c_d * a_p * 0.5*rho_g*(XP(end,:).^2) - QP(end,:,1,1)./QP(end,:,1,2) );
        case 2
            H(2)=plot(YP(end,:),QP(end,:,1,1)./QP(end,:,1,2),'b.');
            v = abs( c_d * a_p * 0.5*rho_g*(YP(end,:).^2) - QP(end,:,1,1)./QP(end,:,1,2) );
        case 3
            H(2)=plot(ZP(end,:),QP(end,:,1,1)./QP(end,:,1,2),'b.');
            v = abs( c_d * a_p * 0.5*rho_g*(ZP(end,:).^2) - QP(end,:,1,1)./QP(end,:,1,2) );
    end

    err = norm(v)/length(v);
    if err>1e-4
        display(['Error: Case ' [ddir,chid{i}] ' error = ' num2str(err)])
    end

end % chid_for

xlabel('Position (m)','FontSize',Label_Font_Size)
ylabel('Drag Force (N)','FontSize',Label_Font_Size)
lh=legend(H,'exact','FDS part');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

Git_Filename = [ddir,chid{1},'_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/part_drag_profile');









