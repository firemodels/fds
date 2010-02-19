% McDermott
% 2-15-10
% radial_profile.m
%
% Plot radial profile from experiment vs. FDS result

function []=radial_profile(exp_file,fds_file,plot_file,devc_col,error_frac,exp_format,fds_format,yaxis_label,xmin,xmax,dx,ymin,ymax,dy)

% experimental data
M = csvread(exp_file,1,0);
r = M(:,1);
w = M(:,2);
e = error_frac*abs(w); % percent error matches DesJardin paper
H(1)=errorbar(r,w,e,exp_format); hold on

% FDS data
M = csvread(fds_file,2,0);
R = xmin:.05:xmax;
W = M(1,devc_col);
H(2)=plot(R,W,fds_format,'LineWidth',2); hold off

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('radial position (m)')
ylabel(yaxis_label)
axis([xmin xmax ymin ymax])
set(gca,'YTick',ymin:dy:ymax)
set(gca,'YMinorTick','on')
set(gca,'XTick',xmin:dx:xmax)
set(gca,'XMinorTick','on')
%set(gca,'YMinorGrid','on')

h = legend(H,'Exp','FDS','Location','Northwest');
set(h,'Interpreter','LaTeX')

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../../Manuals/FDS_5_Validation_Guide/FIGURES/',plot_file])

