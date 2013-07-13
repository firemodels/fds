% McDermott
% 12-28-11
% jet_decay.m

close all
clear all

addpath ../../Validation/Turbulent_Jet/FDS_Output_Files/

% gather FDS results

M = importdata('jet_csmag_dx10cm_line.csv',',',2);     u_csmag_10     = M.data(:,2);
M = importdata('jet_dsmag_dx10cm_line.csv',',',2);     u_dsmag_10     = M.data(:,2);
M = importdata('jet_deardorff_dx10cm_line.csv',',',2); u_deardorff_10 = M.data(:,2);
M = importdata('jet_vreman_dx10cm_line.csv',',',2);    u_vreman_10    = M.data(:,2);

M = importdata('jet_csmag_dx5cm_line.csv',',',2);      u_csmag_5      = M.data(:,2);
M = importdata('jet_dsmag_dx5cm_line.csv',',',2);      u_dsmag_5      = M.data(:,2);
M = importdata('jet_deardorff_dx5cm_line.csv',',',2);  u_deardorff_5  = M.data(:,2);
M = importdata('jet_vreman_dx5cm_line.csv',',',2);     u_vreman_5     = M.data(:,2);

% analytical solutions

x = M.data(:,1);
u_0 = 2.1;
h = 0.8;
b = 0.8;
m_1 = 0.12;
m_2 = 0.20;

for j=1:length(x)
	if x(j)<h/m_1
		u_1(j) = u_0;
	else
        u_1(j) = u_0/(m_1*x(j))*sqrt(b*h);
	end
end

for j=1:length(x)
	if x(j)<h/m_2
		u_2(j) = u_0;
	else
        u_2(j) = u_0/(m_2*x(j))*sqrt(b*h);
	end
end

plot_style

H(1)=plot(x/h,u_1/u_0,'k--'); hold on
H(2)=plot(x/h,u_2/u_0,'k-');
H(3)=plot(x/h,u_csmag_10/u_0,'--','color',[0 .5 0]);
H(4)=plot(x/h,u_csmag_5/u_0,'-','color',[0 .5 0]);
H(5)=plot(x/h,u_dsmag_10/u_0,'c--');
H(6)=plot(x/h,u_dsmag_5/u_0,'c-');
H(7)=plot(x/h,u_deardorff_10/u_0,'b--');
H(8)=plot(x/h,u_deardorff_5/u_0,'b-');
H(9)=plot(x/h,u_vreman_10/u_0,'r--');
H(10)=plot(x/h,u_vreman_5/u_0,'r-');

axis([0 25 0.2 1.2])
set(gca,'XTick',[0 5 10 15 20 25])
set(gca,'YTick',[0.2 0.4 0.6 0.8 1.0 1.2])
Paper_Width=1.4*Paper_Width;
Plot_Width =1.4*Plot_Width;
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

text(1,1.1,'Jet Centerline Velocity Decay','FontSize',Label_Font_Size,'FontName',Font_Name)

xlabel('{\it x/h}','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter);
ylabel('{\it u}_{max}/{\it u}_0','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
legend_handle = legend(H,'analytical, {\it m}=0.12','analytical, {\it m}=0.20',...
	                     'csmag, {\it h/\deltax}=8','csmag, {\it h/\deltax}=16',...
	                     'dsmag, {\it h/\deltax}=8','dsmag, {\it h/\deltax}=16',...
		                 'Deardorff, {\it h/\deltax}=8','Deardorff, {h/\deltax}=16',...
		                 'Vreman, {\it h/\deltax}=8','Vreman, {\it h/\deltax}=16',...
						 'Location','EastOutside');
set(legend_handle,'Interpreter',Font_Interpreter);
set(legend_handle,'Fontsize',Key_Font_Size);
set(legend_handle,'Box','on');

% add SVN if file is available

SVN_Filename = ['jet_dsmag_dx5cm_svn.txt'];
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
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/jet_decay'])


