% McDermott
% 10-30-13
% backward_facing_step.m

close all
clear all

expdir = '../../Validation/Backward_Facing_Step/Experimental_Data/';
datdir = '../../Validation/Backward_Facing_Step/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

plot_style

nx = [5 10 20];
lnx = length(nx);
rho = 1.1992662; % density, kg/m^3
U_0 = 7.72;      % reference velocity, m/s
h = 0.0098;      % step height, m
dx = h./nx;
fds_marker = {'b-' 'r-' 'm-'};
fds_key = {'FDS {\it h/\deltax}=5' 'FDS {\it h/\deltax}=10' 'FDS {\it h/\deltax}=20'};

if ~exist([expdir,'backward_facing_step_data.csv'])
    display(['Error: File ' [datdir,'backward_facing_step_data.csv'] ' does not exist. Skipping case.'])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% streamwise data along bottom of channel for determining reattachment location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot experimental data

j = find(strcmp(D.colheaders,'Cf-x/h'));
xoh = D.data(:,j);
j = find(strcmp(D.colheaders,'Cf'));
Cf = D.data(:,j);

figure
H(1)=plot(xoh,Cf,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'u_tau-x'));
    x = M{i}.data(:,j);
    I = find(x>0);
    x = x(I);
    
    j = find(strcmp(M{i}.colheaders,'u_tau'));
    u_tau = M{i}.data(I,j);
    tau_w = rho*u_tau.^2;

    j = find(strcmp(M{i}.colheaders,'u_wall'));
    u_wall = M{i}.data(I,j);

    Cf_fds = sign(u_wall).*tau_w/(.5*rho*U_0^2);
    
    H(1+i)=plot(x/h,Cf_fds,fds_marker{i}); hold on
end

ylabel('{\it C}_f','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh = xlabel('\it x/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)*1.025    Pos(3)])

%text(.9525,5.55,'Inlet Velocity','Interpreter',Font_Interpreter,'Fontsize',Title_Font_Size,'FontName',Font_Name)

plot([0 20], [0 0], 'k--'); 

set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

lh = legend(H,'Exp data',fds_key{1:lnx},'Location','Southeast');
set(lh,'Interpreter',Font_Interpreter)
legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_Cf'])


% pressure coefficient along bottom of channel

% plot experimental data

j = find(strcmp(D.colheaders,'Cp-x/h'));
xoh = D.data(:,j);
I = find(xoh<=20);
xoh = xoh(I);
j = find(strcmp(D.colheaders,'Cp'));
Cp = D.data(I,j);

figure
H(1)=plot(xoh,Cp,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'cp-x'));
    x = M{i}.data(:,j);
    I = find(x>0);
    x = x(I);
    
    j = find(strcmp(M{i}.colheaders,'cp'));
    Cp_fds = M{i}.data(I,j)/(U_0^2);
    Cp_fds = Cp_fds + Cp(end) - Cp_fds(end);
    
    H(1+i)=plot(x/h,Cp_fds,fds_marker{i}); hold on
end

ylabel('{\it C}_p','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh = xlabel('\it x/h','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)*1.025    Pos(3)])

%text(.9525,5.55,'Inlet Velocity','Interpreter',Font_Interpreter,'Fontsize',Title_Font_Size,'FontName',Font_Name)

set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

lh = legend(H,'Exp data',fds_key{1:lnx},'Location','Southeast');
set(lh,'Interpreter',Font_Interpreter)
legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_Cp'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal streamwise velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,4,1)

j = find(strcmp(D.colheaders,'z -3'));
z = D.data(:,j);
I = find(z>0);
z=z(I)*0.001+h;
j = find(strcmp(D.colheaders,'U -3'));
u_data = D.data(I,j);

H(1)=plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'U-VEL -3-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'U-VEL -3'));
    u_fds = M{i}.data(I,j);
    
    H(1+i)=plot(u_fds,z,fds_marker{i});
end

axis([0 10 0 3.5*h])

ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh=xlabel('{\it U} (m/s)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)*5    Pos(2)    Pos(3)])

subplot(1,4,2)

j = find(strcmp(D.colheaders,'z 4'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'U 4'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'U-VEL 4-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'U-VEL 4'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,3);

j = find(strcmp(D.colheaders,'z 6'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'U 6'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'U-VEL 6-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'U-VEL 6'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i}); hold on
end

axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,4);

j = find(strcmp(D.colheaders,'z 10'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'U 10'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'U-VEL 10-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'U-VEL 10'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

axis([0 10 0 3.5*h])
set(gca,'YTickLabel',[])

% lh = legend(H,'Exp data',fds_key{1:lnx},'Location','SouthEast');
% set(lh,'Interpreter',Font_Interpreter)
% legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_U'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal vertical velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,4,1)

j = find(strcmp(D.colheaders,'z -3'));
z = D.data(:,j);
I = find(z>0);
z=z(I)*0.001+h;
j = find(strcmp(D.colheaders,'W -3'));
u_data = D.data(I,j);

H(1)=plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'W-VEL -3-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'W-VEL -3'));
    u_fds = M{i}.data(I,j);
    
    H(1+i)=plot(u_fds,z,fds_marker{i});
end

axis([-.1 .1 0 3.5*h])

ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh=xlabel('{\it W} (m/s)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)    Pos(3)])

subplot(1,4,2)

j = find(strcmp(D.colheaders,'z 4'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'W 4'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'W-VEL 4-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'W-VEL 4'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,3);

j = find(strcmp(D.colheaders,'z 6'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'W 6'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'W-VEL 6-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'W-VEL 6'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i}); hold on
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,4);

j = find(strcmp(D.colheaders,'z 10'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'W 10'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'W-VEL 10-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'W-VEL 10'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([0 10 0 3.5*h])
set(gca,'YTickLabel',[])

% lh = legend(H,'Exp data',fds_key{1:lnx},'Location','SouthEast');
% set(lh,'Interpreter',Font_Interpreter)
% legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_W'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal uu profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,4,1)

j = find(strcmp(D.colheaders,'z -3'));
z = D.data(:,j);
I = find(z>0);
z=z(I)*0.001+h;
j = find(strcmp(D.colheaders,'uu -3'));
u_data = D.data(I,j);

H(1)=plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uu -3-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uu -3'));
    u_fds = M{i}.data(I,j);
    
    H(1+i)=plot(u_fds,z,fds_marker{i});
end

axis([min(min(u_data),min(u_fds)) max(max(u_data),max(u_fds)) 0 3.5*h])

ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh=xlabel('{\it uu} (m^2/s^2)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)    Pos(3)])

subplot(1,4,2)

j = find(strcmp(D.colheaders,'z 4'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uu 4'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uu 4-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uu 4'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,3);

j = find(strcmp(D.colheaders,'z 6'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uu 6'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uu 6-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uu 6'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i}); hold on
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,4);

j = find(strcmp(D.colheaders,'z 10'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uu 10'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uu 10-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uu 10'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([0 10 0 3.5*h])
set(gca,'YTickLabel',[])

% lh = legend(H,'Exp data',fds_key{1:lnx},'Location','SouthEast');
% set(lh,'Interpreter',Font_Interpreter)
% legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_uu'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal ww profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,4,1)

j = find(strcmp(D.colheaders,'z -3'));
z = D.data(:,j);
I = find(z>0);
z=z(I)*0.001+h;
j = find(strcmp(D.colheaders,'ww -3'));
u_data = D.data(I,j);

H(1)=plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'ww -3-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'ww -3'));
    u_fds = M{i}.data(I,j);
    
    H(1+i)=plot(u_fds,z,fds_marker{i});
end

axis([min(min(u_fds),min(u_data)) max(max(u_fds),max(u_data)) 0 3.5*h])

ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh=xlabel('{\it ww} (m^2/s^2)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)    Pos(3)])

subplot(1,4,2)

j = find(strcmp(D.colheaders,'z 4'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'ww 4'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'ww 4-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'ww 4'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,3);

j = find(strcmp(D.colheaders,'z 6'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'ww 6'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'ww 6-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'ww 6'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i}); hold on
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,4);

j = find(strcmp(D.colheaders,'z 10'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'ww 10'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'ww 10-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'ww 10'));
    u_fds = M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([0 10 0 3.5*h])
set(gca,'YTickLabel',[])

% lh = legend(H,'Exp data',fds_key{1:lnx},'Location','SouthEast');
% set(lh,'Interpreter',Font_Interpreter)
% legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_ww'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall-normal uw profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,4,1)

j = find(strcmp(D.colheaders,'z -3'));
z = D.data(:,j);
I = find(z>0);
z=z(I)*0.001+h;
j = find(strcmp(D.colheaders,'uw -3'));
u_data = D.data(I,j);

H(1)=plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uw -3-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uw -3'));
    u_fds = -M{i}.data(I,j);
    
    H(1+i)=plot(u_fds,z,fds_marker{i});
end

axis([min(min(u_fds),min(u_data)) max(max(u_fds),max(u_data)) 0 3.5*h])

ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name)
xh=xlabel('{\it -uw} (m^2/s^2)','Interpreter',Font_Interpreter,'Fontsize',Label_Font_Size,'FontName',Font_Name);
Pos=get(xh,'Position');
set(xh,'Position',[Pos(1)    Pos(2)    Pos(3)])

subplot(1,4,2)

j = find(strcmp(D.colheaders,'z 4'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uw 4'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uw 4-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uw 4'));
    u_fds = -M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,3);

j = find(strcmp(D.colheaders,'z 6'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uw 6'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uw 6-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uw 6'));
    u_fds = -M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i}); hold on
end

%axis([-5 10 0 3.5*h])
set(gca,'YTickLabel',[])

subplot(1,4,4);

j = find(strcmp(D.colheaders,'z 10'));
z = D.data(:,j)*0.001;
j = find(strcmp(D.colheaders,'uw 10'));
u_data = D.data(:,j);

plot(u_data,z,'ko'); hold on

for i=1:lnx
    j = find(strcmp(M{i}.colheaders,'uw 10-z'));
    z = M{i}.data(:,j);
    I = find(z>0);
    z = z(I);
    
    j = find(strcmp(M{i}.colheaders,'uw 10'));
    u_fds = -M{i}.data(I,j);
    
    plot(u_fds,z,fds_marker{i})
end

%axis([0 10 0 3.5*h])
set(gca,'YTickLabel',[])

% lh = legend(H,'Exp data',fds_key{1:lnx},'Location','SouthEast');
% set(lh,'Interpreter',Font_Interpreter)
% legend boxoff

% add SVN if file is available

SVN_Filename = [datdir,'backward_facing_step_',num2str(nx(1)),'_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'backward_facing_step_uw'])
