% Mueller
% 11-01-2022
% level_set_ellipse.m

% This script calculates expected fire spread behavior with applied
% wind and slope, using an assumption of elliptical spread and the
% Rothermel fire spread model (Bova et al., 2016).
% Expected values are compared to FDS predictions using the level set mode

close all
clear all

out_dir = '../../Verification/WUI/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% figure setup
figure
plot_style
colors=['k','r','b','g'];

% simulation parameters
U = [0,5]; % windspeeds tested (m/s)
R0 = 0.05; % zero-wind, zero-slope spread rate (m/s)
theta_w = 90; % compass direction of wind vector
theta_s = 45; % compass direction of increasing slope (aspect + 180)
slope_angle = [0,30]; % degrees of slope tested
theta = 0:45:315; % device angles

for slope=slope_angle
    
    % initialize figure
    clf
    hold on
    box on    

    for ui=1:length(U)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot FDS spread
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % import arrival times
        CHID = [out_dir,...
            'LS_ellipse_',num2str(U(ui)),'ms_',num2str(slope,'%02d'),'deg'];
        % Read the output data
        if ~exist([out_dir,CHID,'_devc.csv'],'file') 
            disp(['Error: File ',CHID,'_devc.csv does not exist. Skipping case.']);
            return;
        end
        DEVC = importdata([out_dir,[CHID,'_devc.csv']],',', 2);
        % distance from ignition for each point (along sloped surface)
        dist = sqrt(4^2+(4*tand(slope)*cosd(theta-theta_s)).^2);
        % spread rate from arrival time at each point
        r = dist./DEVC.data(end,2:end);
        % plot distance traveled in 10 seconds in cartesian space
        plot(10*r.*sind(theta),10*r.*cosd(theta),'o','color',colors(ui),...
            'displayname',['FDS ',num2str(U(ui)),' m/s']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot theoretical spread
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Midflame wind speed from Andrews, default FDS fuel depth of 1 m
        HT=0.2;
        UMF = U(ui)*1.83/log((20+1.18*HT)/(0.43*HT));
        
        % determine adjustment factors for wind and slope
        [phi_w_x,phi_w_y] = wind_adj(UMF,theta_w);
        [phi_s_x,phi_s_y] = ...
            slope_adj(sind(theta_s)*tand(slope),cosd(theta_s)*tand(slope));
        phi_x = phi_w_x + phi_s_x;
        phi_y = phi_w_y + phi_s_y;
    
        % determine virtual wind vector
        [u_v,v_v] = virtual_wind(phi_s_x,phi_s_y);
        u_v = UMF*sind(theta_w) + u_v;
        v_v = UMF*cosd(theta_w) + v_v;
        UMF_v = sqrt(u_v^2 + v_v^2);
        theta_w_v = 90 - atan2d(v_v,u_v); % compass angle of virtual wind vector
    
        % add slope and wind adjustments to ROS
        R = R0*(1 + sqrt(phi_x^2 + phi_y^2));
        
        % elliptical parameters
        LB = 0.936*exp(0.2566*UMF_v) + 0.461*exp(-0.1548*UMF_v) - 0.397;
        LB = max(1,min(LB,8)); %length to bredth ratio from Anderson and Finney
        HB = (LB + sqrt(LB^2-1))/(LB-sqrt(LB^2-1)); %head to backing ratio
        b = 0.5*(R+R/HB);
        a = b/LB;
        c = b-R/HB;
        % all possible directional unit vectors
        x_s = cosd(0:359);
        y_s = sind(0:359);

        % denominator
        D=(a^2*(x_s.*sind(theta_w_v)+y_s.*cosd(theta_w_v)).^2 +...
            b^2*(x_s.*cosd(theta_w_v)-y_s.*sind(theta_w_v)).^2).^-0.5;
        % spread in x-direction
        R_x = D.*(a^2*cosd(theta_w_v).*(x_s.*sind(theta_w_v)+y_s.*cosd(theta_w_v)) - ...
            b^2*sind(theta_w_v).*(x_s.*cosd(theta_w_v)-y_s.*sind(theta_w_v))) + c*sind(theta_w_v);
        % spread in y-direction
        R_y = D.*(-a^2*sind(theta_w_v).*(x_s.*sind(theta_w_v)+y_s.*cosd(theta_w_v)) - ...
            b^2*cosd(theta_w_v).*(x_s.*cosd(theta_w_v)-y_s.*sind(theta_w_v))) + c*cosd(theta_w_v);
        % distance traveled in 10 seconds
        plot(10*R_x,10*R_y,'-','color',colors(ui),...
            'displayname',['Expected ',num2str(U(ui)),' m/s']);
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute relative error in spread rate for comparison points
        R_target = sqrt(R_x.^2 + R_y.^2);
        theta_target = 90 - atan2d(R_y,R_x);
        theta_target(theta_target<0) =  360 + theta_target(theta_target<0);
        [theta_target,idx] = sort(theta_target,'ascend');
        R_target = R_target(idx);
        R_target = interp1(theta_target,R_target,theta,'linear','extrap');
        % either create or append to table of error data
        if ~exist('error_table','var')
            error_table=table(theta',[abs(r-R_target)./R_target]',...
                'VariableNames',{'angle',...
                [num2str(U(ui)),' m/s; ',num2str(slope,'%02d'),'°']});
        else
            error_table.([num2str(U(ui)),' m/s; ',num2str(slope,'%02d'),'°'])=...
                (abs(r-R_target)./R_target)';
        end

    end
  
    % plot ignition point
    plot(0,0,'r.','HandleVisibility','off');
    % standardize figure
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Scat_Label_Font_Size)    
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    
    % label axes
    axis([-2 10 -4 4])
    xlabel('x (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('y (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)     
    % add legend
    lh=legend();
    set(lh,'FontSize',Key_Font_Size,'location','southeast')
    % add annotations for wind and slope direction
    annotation('textarrow',[0.25,0.35],[0.2,0.2],'String','wind',...
        'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
    if slope>0
        annotation('textarrow',[0.45,0.45+0.1*cosd(theta_s)],[0.23,0.23+0.1*sind(theta_s)],...
            'String','slope','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
    end

   % add Git revision
   Git_Filename = [out_dir,CHID,'_git.txt'];
   addverstr(gca,Git_Filename,'linear')

    % save figure
    print(gcf,Image_File_Type,...
        [plot_dir,'level_set_ellipse_',num2str(slope,'%02d'),'deg'])

end

max_err = max(error_table{:,2:end}(:));
if max_err>.15
    display(['Matlab Warning: LS_ellipse is out of tolerance. Max error = ',num2str(max_err)])
end

% plot error_table to be used for verification
clf
hold on
box on
% standardize figure
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)    
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
for i=2:width(error_table)
    plot(error_table{:,1},error_table{:,i},'.-','color',colors(i-1),...
        'displayname',error_table.Properties.VariableNames{i})
end
xlabel('compass angle (°)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('relative error (-)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)     
lh=legend();

% add Git revision
Git_Filename = [out_dir,CHID,'_git.txt'];
addverstr(gca,Git_Filename,'linear')

% save figure
print(gcf,Image_File_Type,[plot_dir,'level_set_ellipse_error'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fire spread functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate wind adjustment vector
function [phi_w_x,phi_w_y]=wind_adj(U,theta_w)
    U=U*60; %convert to m/min
    sigma=50; % svr (cm^-1), FDS default
    beta=0.01; % packing ratio (-), FDS default
    % Rothermel optimum packing ratio
    beta_op=0.20395*(sigma^-0.8189);
    beta_ratio=beta./beta_op;
    % Rothermel wind coefficient
    C=7.47*exp(-0.8711*sigma.^0.55);
    B=0.15988*sigma.^0.54;
    E=0.715*exp(-0.01094*sigma);
    phi_w_x=C.*(3.281*U).^B.*(beta_ratio).^-E*sind(theta_w);
    phi_w_y=C.*(3.281*U).^B.*(beta_ratio).^-E*cosd(theta_w);
end

% calculate slope adjustment vector
function [phi_s_x,phi_s_y] = slope_adj(dzdx,dzdy)
    beta=0.01; % packing ratio (-), FDS default
    % Rothermel optimum packing ratio
    dzds = sqrt(dzdx^2+dzdy^2);
    phi_s_x = 5.275*beta^-0.3*dzdx*dzds;
    phi_s_y = 5.275*beta^-0.3*dzdy*dzds;
end

% calculate virtual wind vector for slope
function [u_virtual,v_virtual]= virtual_wind(phi_s_x,phi_s_y)
    sigma = 50; % svr (cm^-1), FDS default
    beta = 0.01; % packing ratio (-), FDS default
    % Rothermel optimum packing ratio
    beta_op = 0.20395*(sigma^-0.8189);
    beta_ratio = beta./beta_op;
    % Rothermel wind coefficient
    C = 7.47*exp(-0.8711*sigma.^0.55);
    B = 0.15988*sigma.^0.54;
    E = 0.715*exp(-0.01094*sigma);
    phi_s=sqrt(phi_s_x^2+phi_s_y^2);
    if (phi_s>0)
        uv_tmp = 0.3048/phi_s*(phi_s/C*beta_ratio^E)^(1/B);
    else
        uv_tmp = 0;
    end
    % divide by 60 to maintain units of m/s elsewhere
    u_virtual = uv_tmp*phi_s_x/60;
    v_virtual = uv_tmp*phi_s_y/60;
end
