% McDermott
% 3-5-10
% velocity_signal.m

function []=velocity_signal(devc_file,devc_col,vel_style,tmin,tmax,vmin,vmax,xaxis_title,yaxis_title, ...
                            title_label,text_label,signal_file,svn_file)

plotdir = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/'];
datadir = ['../../Validation/Sandia_Plumes/FDS_Output_Files/'];

M = csvread([datadir,devc_file],2,0);

range = find(M(:,1)>tmin & M(:,1)<tmax);

t = M(range,1);
W = M(range,devc_col);

plot_style

plot(t,W,vel_style)

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel(xaxis_title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(yaxis_title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([tmin tmax vmin vmax ])

xt = tmin + .05*(tmax-tmin);
yt = vmin + .92*(vmax-vmin);
text(xt,yt,title_label,'FontSize',16,'Interpreter',Font_Interpreter)
xt = tmin + .05*(tmax-tmin);
yt = vmin + .84*(vmax-vmin);
text(xt,yt,text_label,'FontSize',16,'Interpreter',Font_Interpreter)

set(gca,'YTick',vmin:1:vmax)
%set(gca,'YMinorTick','on')
set(gca,'XTick',tmin:1:tmax)
%set(gca,'XMinorTick','on')

% add SVN if file is available

if exist(svn_file,'file')
    SVN = importdata(svn_file);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plotdir,'Sandia_Plumes/',signal_file])

