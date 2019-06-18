%!/usr/bin/matlab
%Vanella rotcube_cc_mms_error.m:
%05-24-2019

close all
clear all

% Rotated cube parameters:
mu = 0.01;
D  = 0.01;
rho= 1.00;
L  =   pi; % Cube side length.

nwave=1.;  % Wave number on analytical solution. 
gam  =pi/2;
Az   =0.1;
meanz=0.15;
displ =pi;
dispxy=-1/2*displ*ones(2,1);

% Analytical solutions:
ann_z = @(x,y,t) ...
Az/3.*sin(t)*(1.-cos(2.*nwave*(x-gam)))*(1.-cos(2.*nwave*(y-gam))) - ...
Az/3.*sin(t)+meanz;

ann_u = @(x,y,t) ...
-sin(t)*sin(nwave*x)^2 * sin(2.*nwave*y);

ann_v = @(x,y,t) ...
sin(t)*sin(2.*nwave*x) * sin(nwave*y)^2;

ann_p = @(x,y,t) ...
sin(t)/4.*(2.+cos(2.*nwave*x))*(2.+cos(2.*nwave*y))-sin(t);

ann_H = @(x,y,t) ...
(ann_u(x,y,t)^2+ann_v(x,y,t)^2)/2 + ann_p(x,y,t)/rho;

datadir = '../../Verification/Complex_Geometry/';

ifile_s=1;
ifile_f=3;
ifile_ang = [0; atan(1./2.); atan(1.)];
ifile_str = {'0deg_';'27deg_';'45deg_'};
jfile_str = {'exp';'imp';'obs'};

skip_case = 0;
nfile     = 0;
for ifile=ifile_s:ifile_f

   jfile_s = 1;
   jfile_f = 2;
   if(ifile==1)
       jfile_f = 3;
   end
   for jfile=jfile_s:jfile_f
       nfile = nfile+1;
       file(nfile).name = { ['rotated_cube_' ifile_str{ifile}  '32_' jfile_str{jfile} '_mms.csv'], ...
                            ['rotated_cube_' ifile_str{ifile}  '64_' jfile_str{jfile} '_mms.csv'], ...
                            ['rotated_cube_' ifile_str{ifile} '128_' jfile_str{jfile} '_mms.csv'], ...
                            ['rotated_cube_' ifile_str{ifile} '256_' jfile_str{jfile} '_mms.csv'], ...};
                            ['rotated_cube_' ifile_str{ifile} '384_' jfile_str{jfile} '_mms.csv']};
       file(nfile).nameout = ['rotated_cube_' ifile_str{ifile} jfile_str{jfile} '_mms_connvergence'];
       file(nfile).rotang=ifile_ang(ifile);
       for n=1:length(file(nfile).name)
           if ~exist([datadir,file(nfile).name{n}])
               display(['Error: File ' [datadir,file(nfile).name{n}] ...
                        ' does not exist. Skipping case.'])
               skip_case = 1;
           end
       end
   end
end
if skip_case
    return
end

% Compute errors from *_mms.csv files:
for ifile=1:nfile
    
    e_z = zeros(length(file(ifile).name),1);
    e_u = zeros(length(file(ifile).name),1);
    e_v = zeros(length(file(ifile).name),1);
    e_H = zeros(length(file(ifile).name),1);
    
    rotangle   = file(ifile).rotang;
    ROTMAT     = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
    TROTMAT    = ROTMAT';
    
    for n=1:length(file(ifile).name)
        
%         disp(file(ifile).name{n})
        
        M = importdata([datadir,file(ifile).name{n}],',',1); % FDS results.
        vec = str2num(M.textdata{1,1}); % exact time value
        
        ntot_u = vec(1);
        ntot_v = vec(2);
        ntot_c = vec(3);
        T      = vec(4);
        dx     = vec(5);
        dy     = vec(6);
        
        file(ifile).dx(n) = dx;
        
        % inialize error arrays
        e_z_vec = zeros(ntot_c,1);        
        e_u_vec = zeros(ntot_u,1);
        e_v_vec = zeros(ntot_v,1);
        e_H_vec = zeros(ntot_c,1);
        
        
        p = 0;
        % First u:
        Atot_u = 0;
        for i=1:ntot_u
            p=p+1;
            iscf = M.data(p,1);
            xglob= M.data(p,2);
            yglob= M.data(p,3);
            area = M.data(p,4);
            uglob= M.data(p,5);
            
            % Compute analytical value of u in global coordinate system:
            XGLOB = [xglob yglob]' - displ;
            XLOC  = TROTMAT*XGLOB  - dispxy;
            
            % Analytical Velocity in local axes:
            ULOC(1,1) = ann_u(XLOC(1),XLOC(2),T);
            ULOC(2,1) = ann_v(XLOC(1),XLOC(2),T);
            
            % Analytical Velocity in global axes:
            UGLOB = ROTMAT*ULOC;
            
            uglob_ann = UGLOB(1);
            
            % Error:
            e_u_vec(i) = uglob - uglob_ann;
            
            Atot_u = Atot_u + area;
        end
        
        file(ifile).e_u_1(n)   = norm(e_u_vec(1:ntot_u),1)/ntot_u;
        file(ifile).e_u_2(n)   = norm(e_u_vec(1:ntot_u),2)/sqrt(ntot_u);
        file(ifile).e_u_i(n)   = norm(e_u_vec(1:ntot_u),inf);
        
        % Second v:
        Atot_v = 0;
        for i=1:ntot_v
            p=p+1;
            iscf = M.data(p,1);
            xglob= M.data(p,2);
            yglob= M.data(p,3);
            area = M.data(p,4);
            vglob= M.data(p,5);
            
            % Compute analytical value of v in global coordinate system:
            XGLOB = [xglob yglob]' - displ;
            XLOC  = TROTMAT*XGLOB  - dispxy;
            
            % Analytical Velocity in local axes:
            ULOC(1,1) = ann_u(XLOC(1),XLOC(2),T);
            ULOC(2,1) = ann_v(XLOC(1),XLOC(2),T);
            
            % Analytical Velocity in global axes:
            UGLOB = ROTMAT*ULOC;
            
            vglob_ann = UGLOB(2);
            
            % Error:
            e_v_vec(i) = vglob - vglob_ann;
            
            Atot_v = Atot_v + area;
        end
        
        file(ifile).e_v_1(n)   = norm(e_v_vec(1:ntot_v),1)/ntot_v;
        file(ifile).e_v_2(n)   = norm(e_v_vec(1:ntot_v),2)/sqrt(ntot_v);
        file(ifile).e_v_i(n)   = norm(e_v_vec(1:ntot_v),inf);

        % Finally cell centered variables:
        Vtot_c = 0;
        DELP   = 0;
        p_c    = p;
        for i=1:ntot_c
            p=p+1;
            iscf = M.data(p,1);
            xglob= M.data(p,2);
            yglob= M.data(p,3);
            vol  = M.data(p,4);        
            z    = M.data(p,5); 
            H    = M.data(p,6);
            pres = M.data(p,7);
            
            % Compute analytical scalars in global coordinate system:
            XGLOB = [xglob yglob]' - displ;
            XLOC  = TROTMAT*XGLOB  - dispxy;
            
            % Analytical Velocity in local axes:
            z_ann = ann_z(XLOC(1),XLOC(2),T);
            H_ann = ann_H(XLOC(1),XLOC(2),T);
            p_ann = ann_p(XLOC(1),XLOC(2),T);
            
            if(i==1) 
                DELP = pres - p_ann;
            end
            
            % Error:
            e_z_vec(i) = z - z_ann;
            e_H_vec(i) = H - H_ann - DELP/rho;
            e_p_vec(i) = pres - p_ann - DELP;
            
            Vtot_c = Vtot_c + vol;
        end
          
        file(ifile).e_z_1(n)   = norm(e_z_vec(1:ntot_c),1)/ntot_c;
        file(ifile).e_z_2(n)   = norm(e_z_vec(1:ntot_c),2)/sqrt(ntot_c);
        file(ifile).e_z_i(n)   = norm(e_z_vec(1:ntot_c),inf);

        file(ifile).e_H_1(n)   = norm(e_H_vec(1:ntot_c),1)/ntot_c;
        file(ifile).e_H_2(n)   = norm(e_H_vec(1:ntot_c),2)/sqrt(ntot_c);
        file(ifile).e_H_i(n)   = norm(e_H_vec(1:ntot_c),inf);
       
        file(ifile).e_p_1(n)   = norm(e_p_vec(1:ntot_c),1)/ntot_c;
        file(ifile).e_p_2(n)   = norm(e_p_vec(1:ntot_c),2)/sqrt(ntot_c);
        file(ifile).e_p_i(n)   = norm(e_p_vec(1:ntot_c),inf);
        [dummy,i] = max(abs(e_p_vec(1:ntot_c)));
        file(ifile).e_p_i_loc(n) = i;
        file(ifile).e_p_i_xy(n,1:3) = [ M.data(p,7) M.data(p_c+i,2) M.data(p_c+i,3) ];

    end
    
    % Error orders:
%     p_L2_u = log(file(ifile).e_u_2(1:end-1)./file(ifile).e_u_2(2:end))./...
%              log(file(ifile).dx(1:end-1)  ./file(ifile).dx(2:end));
%     p_Li_u = log(file(ifile).e_u_i(1:end-1)./file(ifile).e_u_i(2:end))./...
%              log(file(ifile).dx(1:end-1)  ./file(ifile).dx(2:end));
%     disp(' ')
%     disp(['L2 Err u p:' file(ifile).nameout])
%     disp('    dx        L2e        p2          Linfe      pinf')
%     for ip=1:length(p_L2_u)
%        disp([num2str(file(ifile).dx(ip+1))    ', ' ...
%              num2str(file(ifile).e_u_2(ip+1)) ', ' ...
%              num2str(p_L2_u(ip))              ', ' ...
%              num2str(file(ifile).e_u_i(ip+1)) ', ' ...
%              num2str(p_Li_u(ip))])
%     end
%     disp(' ')
%     p_L2_z = log(file(ifile).e_z_2(1:end-1)./file(ifile).e_z_2(2:end))./...
%              log(file(ifile).dx(1:end-1)  ./file(ifile).dx(2:end));
%     p_Li_z = log(file(ifile).e_z_i(1:end-1)./file(ifile).e_z_i(2:end))./...
%              log(file(ifile).dx(1:end-1)  ./file(ifile).dx(2:end));
%     disp(' ')
%     disp(['L2 Err z p:' file(ifile).nameout])
%     disp('    dx        L2e        p2          Linfe      pinf')
%     for ip=1:length(p_L2_z)
%        disp([num2str(file(ifile).dx(ip+1))    ', ' ...
%              num2str(file(ifile).e_z_2(ip+1)) ', ' ...
%              num2str(p_L2_z(ip))              ', ' ...
%              num2str(file(ifile).e_z_i(ip+1)) ', ' ...
%              num2str(p_Li_z(ip))])
%     end
%     disp(' ')

    plot_style
    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    
    hh(1)=loglog(file(ifile).dx,file(ifile).e_z_2,'ks-');
    hold on
    hh(2)=loglog(file(ifile).dx,file(ifile).e_u_2,'k>-');
    hh(3)=loglog(file(ifile).dx,file(ifile).e_H_2,'k+-');
    hh(4)=loglog(file(ifile).dx,10^-1*file(ifile).dx,'k--');
    hh(5)=loglog(file(ifile).dx,2*10^-3*file(ifile).dx.^2,'k-');
    if(ifile <= 3)
       axis([2*10^-3 3*10^-1 10^-7 10^-1]) 
    elseif(ifile <=5)
       axis([10^-2 5*10^-1 10^-7 7*10^-1])
    else
       axis([5*10^-3 5*10^-1 10^-7 10^-1])
    end
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Title_Font_Size)
    
    xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',...
           Font_Interpreter,'Fontname',Font_Name)
    ylabel('{\it L_2} Error','FontSize',Title_Font_Size,'Interpreter', ...
           Font_Interpreter,'Fontname',Font_Name)
    lh=legend(hh,'FDS {\it Z}','FDS {\it u}',...
              'FDS {\it H}','{\it O(\Deltax)}','{\it O(\Deltax^2)}',...
              'location','northwest');
    set(lh,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,...
        'Fontname',Font_Name)
    
    % add Git version if file is available
    Git_Filename = [datadir,'rotated_cube_45deg_256_exp_git.txt'];
    addverstr(gca,Git_Filename,'loglog')
    
    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    
    strng=['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/' ...
           file(ifile).nameout];
    
    print(gcf,'-dpdf',strng)
        
    if (ifile==3) % OBST case:
        % check errors
        if file(ifile).e_z_2(end) > 2e-6
            display(['Matlab Warning: Neumann BC Species in rotated_cube OBST ' ...
            file(ifile).name{end} ' is out of tolerance. e_z = ',      ...
            file(ifile).e_z_2(end)])
        end
        if file(ifile).e_u_2(end) > 2e-5
            display(['Matlab Warning: Velocity in rotated_cube OBST '             ...
                     file(ifile).name{end} ' is out of tolerance. e_u = ', ...
                     num2str(file(ifile).e_u_2(end))])
        end
        if file(ifile).e_H_2(end) > 1e-3
            display(['Matlab Warning: Pressure in rotated_cube OBST '             ...
                     file(ifile).name{end} ' is out of tolerance. e_H = ', ...
                     num2str(file(ifile).e_H_2(end))])
        end
    elseif(ifile==4)
        % check errors
        if file(ifile).e_z_2(end) > 6e-5
            display(['Matlab Warning: Neumann BC Species in rotated_cube 27deg exp ' ...
            file(ifile).name{end} ' is out of tolerance. e_z = ',      ...
            file(ifile).e_z_2(end)])
        end
        if file(ifile).e_u_2(end) > 1.5e-4
            display(['Matlab Warning: Velocity in rotated_cube 27deg exp '             ...
                     file(ifile).name{end} ' is out of tolerance. e_u = ', ...
                     num2str(file(ifile).e_u_2(end))])
        end
        if file(ifile).e_H_2(end) > 2.5e-3
            display(['Matlab Warning: Pressure in rotated_cube 27deg exp '             ...
                     file(ifile).name{end} ' is out of tolerance. e_H = ', ...
                     num2str(file(ifile).e_H_2(end))])
        end
    elseif(ifile==5)
        % check errors
        if file(ifile).e_z_2(end) > 2e-5
            display(['Matlab Warning: Neumann BC Species in rotated_cube 27deg imp ' ...
            file(ifile).name{end} ' is out of tolerance. e_z = ',      ...
            file(ifile).e_z_2(end)])
        end
        if file(ifile).e_u_2(end) > 1.5e-4
            display(['Matlab Warning: Velocity in rotated_cube 27deg imp '             ...
                     file(ifile).name{end} ' is out of tolerance. e_u = ', ...
                     num2str(file(ifile).e_u_2(end))])
        end
        if file(ifile).e_H_2(end) > 2.5e-3
            display(['Matlab Warning: Pressure in rotated_cube 27deg imp '             ...
                     file(ifile).name{end} ' is out of tolerance. e_H = ', ...
                     num2str(file(ifile).e_H_2(end))])
        end        
    end
end
