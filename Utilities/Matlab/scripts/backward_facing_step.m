% Toms
% 8-8-14
% backward_facing_step.m

function main()

clear all
close all

disp('backward_facing_step ...')

expdir = '../../../exp/Backward_Facing_Step/';
datdir = '../../../out/Backward_Facing_Step/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Backward_Facing_Step/';

rkappa = 1/.41;
B = 5.2;
rho = 1.1992662; % density, kg/m^3
U_0 = 7.72;      % reference velocity, m/s
h = 0.0098;      % step height, m
dx = h/20;
visc = 1.7801E-5;
visc_nu = visc/rho;
x_loc = {'-3','4','6','10'};

nx = [5 10 20];
lnx = length(nx);
dx = h./nx;
fds_marker = {'b-' 'r-' 'm-'};
fds_key = {'{\it h/\deltax}=5' '{\it h/\deltax}=10' '{\it h/\deltax}=20'};

if ~exist([expdir,'backward_facing_step_data.csv'])
    display(['Error: File ' [expdir,'backward_facing_step_data.csv'] ' does not exist. Skipping case.'])
    return
end
for i=1:lnx
    if ~exist([datdir,'backward_facing_step_',num2str(nx(i)),'_line.csv'])
        display(['Error: File ' [datdir,'backward_facing_step_',num2str(nx(i)),'_line.csv'] ' does not exist. Skipping case.'])
        return
    end
end

% read experimental and FDS data files

D = importdata([expdir,'backward_facing_step_data.csv'],',',1);
for i=1:lnx
    M{i} = importdata([datdir,'backward_facing_step_',num2str(nx(i)),'_line.csv'],',',2);
end

linetype_vect_leg = [{'-'},{'-'},{'-'}];
linetype_vect_noleg = [{'-'},{'-'},{'-'}];
symbol_vect_noleg = [{'p'},{'s'},{'d'}];

line_color_vect={[20/255,168/255,113/225],[186/255,62/255,62/255],[64/255,97/255,191/255]};

plot_style

f1 = figure(1);
set(f1,'Visible',Figure_Visibility)
p1 = tight_subplot(1,1, [.01 .01],[.13 .08],[.11 .019]);

f2 = figure(2);
set(f2,'Visible',Figure_Visibility);
p2 = tight_subplot(1,1, [.01 .01],[.13 .08],[.141 .019]);

f3 = figure(3);
set(f3,'Visible',Figure_Visibility);
sp1 = tight_subplot(1,4, [.01 .01],[.142 .055],[.108 .01]);

f4 = figure(4);
set(f4,'Visible',Figure_Visibility);
sp2 = tight_subplot(1,4, [.01 .01],[.142 .055],[.108 .01]);

f5 = figure(5);
set(f5,'Visible',Figure_Visibility);
sp3 = tight_subplot(1,4, [.01 .01],[.142 .055],[.108 .01]);

f6 = figure(6);
set(f6,'Visible',Figure_Visibility);
sp4 = tight_subplot(1,4, [.01 .01],[.142 .055],[.108 .01]);

f7 = figure(7);
set(f7,'Visible',Figure_Visibility);
sp5 = tight_subplot(1,4, [.01 .01],[.142 .055],[.108 .01]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% streamwise data along bottom of channel for determining reattachment location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f1)

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

%%%Exp Data%%%

j = find(strcmp(D.colheaders,'Cf-x/h'));
xoh = D.data(:,j);
j = find(strcmp(D.colheaders,'Cf'));
Cf = D.data(:,j);

error = .0005 * ones(length(Cf),1);

h_dat=errorbar(xoh,Cf,error,'ko','markersize',10);

hold on

axis([0 20 -4E-3 5E-3])
set(gca,'YTick',[-4E-3:1E-3:5E-3])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

%%%FDS Data%%%

for q = 1:lnx
    j= find(strcmp(M{q}.colheaders, 'u_tau-x'));
    x = M{q}.data(:,j);
    I = find(x > 0);
    x = x(I);

    j = find(strcmp(M{q}.colheaders, 'u_tau'));
    u_tau = M{q}.data(I,j);
    tau_w = rho*u_tau.^2;

    j = find(strcmp(M{q}.colheaders, 'u_wall'));
    u_wall = M{q}.data(I,j);

    Cf_fds = visc*u_wall/(.5*h/nx(q))/(.5*rho*U_0^2);

    plot(x/h,Cf_fds, char(linetype_vect_leg(q)), 'Color',line_color_vect{q}, 'LineWidth', 1);

    h_leg(q) = plot(x(1:round(length(x)/15):end)/h, Cf_fds(1:round(length(Cf_fds)/15):end), char(symbol_vect_noleg(q)), 'Color', line_color_vect{q}, 'LineWidth', 1, 'MarkerSize', 10);

    yh = ylabel('{\it C_f}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
    xh = xlabel('{\it x/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);

    % if q == 1
    %     set(xh,'Units','Pixels')
    %     xh_pos = get(xh,'Position');
    %     xh_pos(1) = xh_pos(1)*1.073;
    %     xh_pos(2) = xh_pos(2)*1.17;
    %     set(xh,'Position',xh_pos)
    % end

    if q == 1
        Pos=get(yh,'Position');
        set(yh,'Position',[Pos(1)*.95 Pos(2)*.5 Pos(3)])
    end
end

axis([0 20 -4E-3 4E-3])

lh = legend([h_dat,h_leg([1:length(h_leg)])], ['J&D',fds_key([1:length(h_leg)])], 'Location', 'SouthEast');
set(lh,'box','off')
%set(lh,'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);

% add Git revision if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_Cf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pressure coefficient along bottom of channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f2)

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

%%%Exp Data%%%
j = find(strcmp(D.colheaders,'Cp-x/h'));
xoh = D.data(:,j);
I = find(xoh<=20);
xoh = xoh(I);
j = find(strcmp(D.colheaders,'Cp'));
Cp = D.data(I,j);

h_dat=plot(xoh,Cp,'ko','markersize',10);

hold on

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

%%%FDS Data%%%

for q = 1:lnx
    j = find(strcmp(M{q}.colheaders, 'cp-x'));
    x = M{q}.data(:,j);
    I = find(x>0);
    x = x(I);

    j = find(strcmp(M{q}.colheaders, 'cp'));
    Cp_fds = M{q}.data(I,j)/(U_0^2);
    Cp_fds = Cp_fds + Cp(end) - Cp_fds(end);

    plot(x/h,Cp_fds, char(linetype_vect_leg(q)), 'Color',line_color_vect{q}, 'LineWidth', 1);

    h_leg2(q) = plot(x(1:round(length(x)/15):end)/h, Cp_fds(1:round(length(Cp_fds)/15):end), char(symbol_vect_noleg(q)), 'Color', line_color_vect{q}, 'LineWidth', 1, 'MarkerSize', 10);

    yh = ylabel('{\it C_p}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
    xh = xlabel('\it x/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);

    Pos=get(yh, 'Position');
    set(yh, 'Position',[Pos(1)*.93 Pos(2) Pos(3)]);
end

lh = legend([h_dat,h_leg2([1:length(h_leg2)])], ['J&D',fds_key([1:length(h_leg2)])], 'Location', 'SouthEast');

% add version string if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_Cp'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal streamwise velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f3);

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

for i = 1:4

    %%%Exp Data%%%

    j = find(strcmp(D.colheaders,strcat({'z '},x_loc{i})));
    z = D.data(:,j);
    I = find(z>0);
    if i == 1
        z=z(I)*0.001+h;
    else
        z=z(I)*0.001;
    end
    j = find(strcmp(D.colheaders,strcat({'U '},x_loc{i})));
    u_data = D.data(I,j);

    set(f3,'CurrentAxes',sp1(i))
    H(1)=plot(u_data/U_0,z/h,'ko','markersize',10);
    hold on

    %%%FDS data%%%

    for q=1:lnx
        j = find(strcmp(M{q}.colheaders, strcat({'U-VEL '},x_loc{i},'-z')));
        z = M{q}.data(:,j);
        I = find(z>0); %ensures no data for z<0 is incorporated...
        z = z(I);

        j = find(strcmp(M{q}.colheaders, strcat({'U-VEL '},x_loc{i})));
        u_vel = M{q}.data(I,j);

        plot(u_vel/U_0,z/h,char(linetype_vect_leg(q)), 'Color',line_color_vect{q}, 'LineWidth', 1)
        plot(u_vel(1:round(length(u_vel)/15):end)/U_0,z(1:round(length(z)/15):end)/h,char(symbol_vect_noleg(q)), 'MarkerSize', 10, 'Color',line_color_vect{q}, 'LineWidth', 1)

        axis([-.2 1.0 0 3.5])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        Pos = get(gca,'Position');
        set(gca,'Position',[Pos(1) Pos(2)+.015 Pos(3) Pos(4)*0.98])

        if i == 1
            ylabel('{\it z/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            xh=xlabel('{\it <u>/U_0}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
        else
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end
        if q == 1
            tpos = [0.01 3.2]; % from get(th,'Position') where th is a title handle
            th = text(tpos(1),tpos(2),strcat('{\itx/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            set(th,'Position',tpos)
        end
    end

end

% add Git revision if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear',0,1.05)

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_U'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal vertical velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f4);

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

for i = 1:4

    j = find(strcmp(D.colheaders,strcat({'z '},x_loc{i})));
    z = D.data(:,j);
    I = find(z>0);
    if i == 1
        z=z(I)*0.001+h;
    else
        z=z(I)*0.001;
    end
    j = find(strcmp(D.colheaders,strcat({'W '},x_loc{i})));
    w_data = D.data(I,j);

    set(f4,'CurrentAxes',sp2(i))
    plot(w_data/U_0,z/h,'ko','markersize',10);
    hold on

    for q=1:lnx
        j = find(strcmp(M{q}.colheaders,strcat({'W-VEL '},x_loc{i},'-z')));
        z = M{q}.data(:,j);
        I = find(z>0);
        z = z(I);

        j = find(strcmp(M{q}.colheaders,strcat({'W-VEL '},x_loc{i})));
        w_vel = M{q}.data(I,j);

        plot(w_vel/U_0,z/h,char(linetype_vect_noleg(q)), 'MarkerSize', 1.5,'Color',line_color_vect{q}, 'LineWidth', 1);
        plot(w_vel(1:round(length(w_vel)/15):end)/U_0,z(1:round(length(z)/15):end)/h, char(symbol_vect_noleg(q)), 'MarkerSize', 10, 'Color',line_color_vect{q}, 'LineWidth', 1)

        axis([-.1 .1 0 3.5])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        Pos = get(gca,'Position');
        set(gca,'Position',[Pos(1) Pos(2)+.015 Pos(3) Pos(4)*0.98])

        if i == 1
            ylabel('{\it z/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            xh=xlabel('{\it <w>/U_0}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
        else
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end

        if q == 1
            tpos = [0.01 3.2]; % from get(th,'Position') where th is a title handle
            th = text(tpos(1),tpos(2),strcat('{\itx/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            set(th,'Position',tpos)
        end
    end

    %th = title(strcat('{\it x/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name, 'HorizontalAlignment', 'Left');

end

% add Git revision if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear',0,1.05)

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_W'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal uu profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f5);

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

for i = 1:4

    j = find(strcmp(D.colheaders,strcat({'z '}, x_loc{i})));
    z = D.data(:,j);
    I = find(z>0);
    if i == 1
        z=z(I)*0.001+h;
    else
        z=z(I)*0.001;
    end
    j = find(strcmp(D.colheaders,strcat({'uu '},x_loc{i})));
    uu_data = D.data(I,j);

    error = uu_data * .15 / U_0^2;

    set(f5,'CurrentAxes',sp3(i))
    plot(uu_data/U_0^2,z/h,'ko','markersize',10);
    hold on
    %herrorbar(uu_data/U_0^2,z/h, error,'ko');

    for q=1:lnx

        j = find(strcmp(M{q}.colheaders,strcat({'uu '},x_loc{i},'-z')));
        z = M{q}.data(:,j);
        I = find(z>0);
        z = z(I);

        j = find(strcmp(M{q}.colheaders,strcat({'uu '},x_loc{i})));
        uu_fds = M{q}.data(I,j);

        plot(uu_fds/U_0^2,z/h, char(linetype_vect_noleg(q)), 'MarkerSize', 2, 'Color',line_color_vect{q}, 'LineWidth', 1);
        plot(uu_fds(1:round(length(uu_fds)/15):end)/U_0^2,z(1:round(length(z)/15):end)/h,char(symbol_vect_noleg(q)), 'MarkerSize', 10, 'Color',line_color_vect{q}, 'LineWidth', 1)

        axis([0 0.04 0 3.5])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        Pos = get(gca,'Position');
        set(gca,'Position',[Pos(1) Pos(2)+.015 Pos(3) Pos(4)*0.98])

        if i == 1
            ylabel('{\it z/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            xh=xlabel('{\it <uu>/U_0^2}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
        else
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
        end
        if q == 1
            tpos = [0.01 3.2]; % from get(th,'Position') where th is a title handle
            th = text(tpos(1),tpos(2),strcat('{\itx/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            set(th,'Position',tpos)
        end
    end

end

% add Git revision if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear',0,1.05)

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_uu'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal ww profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f6);

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

for i = 1:4

    j = find(strcmp(D.colheaders,strcat({'z '},x_loc{i})));
    z = D.data(:,j);
    I = find(z>0);
    if i == 1
        z=z(I)*0.001+h;
    else
        z=z(I)*0.001;
    end
    j = find(strcmp(D.colheaders,strcat({'ww '},x_loc{i})));
    ww_data = D.data(I,j);

    error = ww_data * .15 / U_0^2;

    set(f6,'CurrentAxes',sp4(i))
    plot(ww_data/U_0^2,z/h,'ko','markersize',10);
    hold on
    %herrorbar(ww_data/U_0^2,z/h,error,'ko');

    for q = 1:lnx
        j = find(strcmp(M{q}.colheaders,strcat({'ww '},x_loc{i},'-z')));
        z = M{q}.data(:,j);
        I = find(z>0);
        z = z(I);


        j = find(strcmp(M{q}.colheaders,strcat({'ww '},x_loc{i})));
        ww_fds = M{q}.data(I,j);


        plot(ww_fds/U_0^2,z/h,char(linetype_vect_noleg(q)), 'MarkerSize', 1.5,'Color',line_color_vect{q}, 'LineWidth', 1);
        plot(ww_fds(1:round(length(ww_fds)/15):end)/U_0^2,z(1:round(length(z)/15):end)/h,char(symbol_vect_noleg(q)), 'MarkerSize', 10, 'Color',line_color_vect{q}, 'LineWidth', 1)

        axis([0 .04 0 3.5])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        Pos = get(gca,'Position');
        set(gca,'Position',[Pos(1) Pos(2)+.015 Pos(3) Pos(4)*0.98])

        if i == 1
            ylabel('{\it z/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
            xh=xlabel('{\it <ww>/U_0^2}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
        else
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
        end

        if q == 1
            tpos = [0.01 3.2]; % from get(th,'Position') where th is a title handle
            th = text(tpos(1),tpos(2),strcat('{\itx/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            set(th,'Position',tpos)
        end
    end

end

% add Git if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear',0,1.05)

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_ww'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal uw profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'CurrentFigure', f7);

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

for i = 1:4
    j = find(strcmp(D.colheaders,strcat({'z '},x_loc{i})));
    z = D.data(:,j);
    I = find(z>0);
    if i == 1
        z=z(I)*0.001+h;
    else
        z=z(I)*0.001;
    end
    j = find(strcmp(D.colheaders,strcat({'uw '},x_loc{i})));
    uw_data = D.data(I,j);

    error = uw_data * .15 / U_0^2;

    set(f7,'CurrentAxes',sp5(i))
    h_dat=plot(uw_data/U_0^2,z/h,'ko','markersize',10);
    hold on
    %herrorbar(uw_data/U_0^2,z/h,error,'ko');

    for q = 1:lnx

        j = find(strcmp(M{q}.colheaders,strcat({'uw '},x_loc{i},'-z')));
        z = M{q}.data(:,j);
        I = find(z>0);
        z = z(I);

        j = find(strcmp(M{q}.colheaders,strcat({'uw '},x_loc{i})));
        uw_fds = -1*M{q}.data(I,j);


        plot(uw_fds/U_0^2,z/h,char(linetype_vect_noleg(q)), 'MarkerSize', 1.5,'Color',line_color_vect{q}, 'LineWidth', 1);
        h_leg(q)=plot(uw_fds(1:round(length(uw_fds)/15):end)/U_0^2,z(1:round(length(z)/15):end)/h,char(symbol_vect_noleg(q)), 'MarkerSize', 10, 'Color',line_color_vect{q}, 'LineWidth', 1);

        axis([0 .04 0 3.5])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Label_Font_Size)

        Pos = get(gca,'Position');
        set(gca,'Position',[Pos(1) Pos(2)+.015 Pos(3) Pos(4)*0.98])

        if i == 1
            ylabel('{\it z/h}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
            xh=xlabel('{\it -<uw>/U_0^2}','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
        else
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
        end

        if q == 1
            tpos = [0.01 3.2]; % from get(th,'Position') where th is a title handle
            th = text(tpos(1),tpos(2),strcat('{\itx/h}','=',x_loc{i}),'Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
            set(th,'Position',tpos)
        end
    end

end

% add version string if file is available

Git_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_git.txt'];
addverstr(gca,Git_Filename,'linear',0,1.05)

% print to pdf

print(gcf,'-dpdf',[pltdir,'backward_facing_step_uw'])

end % main()


%THE FOLLOWING CODE IS ALL EXTERNAL SCRIPTS IMBEDDED IN THIS SCRIPT TO
%LIMIT THE TOMS INFLUENCE ON THE REPOSITORY

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh;

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);

    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end

end


function hh = herrorbar(x, y, l, u, symbol)
%HERRORBAR Horizontal Error bar plot.
%   HERRORBAR(X,Y,L,R) plots the graph of vector X vs. vector Y with
%   horizontal error bars specified by the vectors L and R. L and R contain the
%   left and right error ranges for each point in X. Each error bar
%   is L(i) + R(i) long and is drawn a distance of L(i) to the right and R(i)
%   to the right the points in (X,Y). The vectors X,Y,L and R must all be
%   the same length. If X,Y,L and R are matrices then each column
%   produces a separate line.
%
%   HERRORBAR(X,Y,E) or HERRORBAR(Y,E) plots X with error bars [X-E X+E].
%   HERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'. See PLOT for possibilities.
%
%   H = HERRORBAR(...) returns a vector of line handles.
%
%   Example:
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      herrorbar(x,y,e)
%   draws symmetric horizontal error bars of unit standard deviation.
%
%   This code is based on ERRORBAR provided in MATLAB.
%
%   See also ERRORBAR

%   Jos van der Geest
%   email: jos@jasen.nl
%
%   File history:
%   August 2006 (Jos): I have taken back ownership. I like to thank Greg Aloe from
%   The MathWorks who originally introduced this piece of code to the
%   Matlab File Exchange.
%   September 2003 (Greg Aloe): This code was originally provided by Jos
%   from the newsgroup comp.soft-sys.matlab:
%   http://newsreader.mathworks.com/WebX?50@118.fdnxaJz9btF^1@.eea3ff9
%   After unsuccessfully attempting to contact the orignal author, I
%   decided to take ownership so that others could benefit from finding it
%   on the MATLAB Central File Exchange.

if min(size(x))==1,
    npt = length(x);
    x = x(:);
    y = y(:);
    if nargin > 2,
        if ~isstr(l),
            l = l(:);
        end
        if nargin > 3
            if ~isstr(u)
                u = u(:);
            end
        end
    end
else
    [npt,n] = size(x);
end

if nargin == 3
    if ~isstr(l)
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        [m,n] = size(y);
        x(:) = (1:npt)'*ones(1,n);;
    end
end

if nargin == 4
    if isstr(u),
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end

if nargin == 2
    l = y;
    u = y;
    y = x;
    [m,n] = size(y);
    x(:) = (1:npt)'*ones(1,n);;
    symbol = '-';
end

u = abs(u);
l = abs(l);

if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
    error('The sizes of X, Y, L and U must be the same.');
end

tee = (max(y(:))-min(y(:)))/100; % make tee .02 x-distance for error bars
% changed from errorbar.m
xl = x - l;
xr = x + u;
ytop = y + tee;
ybot = y - tee;
n = size(y,2);
% end change

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
% changed from errorbar.m
xb = zeros(npt*9,n);
xb(1:9:end,:) = xl;
xb(2:9:end,:) = xl;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xr;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = y;
yb(5:9:end,:) = y;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ytop;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;
% end change


[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol); hold on
h = [h;plot(x,y,symbol)];

if ~hold_state, hold off; end

if nargout>0, hh = h; end

end

function cell2csv(filename,cellArray,delimiter)
% Writes cell array content into a *.csv file.
%
% CELL2CSV(filename,cellArray,delimiter)
%
% filename      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% delimiter = seperating sign, normally:',' (default)
%
% by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
if nargin<3
    delimiter = ',';
end

datei = fopen(filename,'w');
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)

        var = eval(['cellArray{z,s}']);

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end

        fprintf(datei,var);

        if s ~= size(cellArray,2)
            fprintf(datei,[delimiter]);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);
end
