% McDermott
% 2-15-10
% radial_profile.m
%
% Plot radial profile from experiment vs. FDS results (up to 4)
%
% Example:
%
% function []=radial_profile(plot_file,devc_col,error_frac, ...
%                            xaxis_label,yaxis_label, ...
%                            rmin,rmax,dr,xmin,xmax,dx,ymin,ymax,dy, ...
%                            exp_file,exp_format,exp_label, ...
%                            fds_file1,fds_format1,fds_label1, ...
%                            fds_file2,fds_format2,fds_label2, ...
%                            fds_file3,fds_format3,fds_label3, ...
%                            fds_file4,fds_format4,fds_label4)

function []=radial_profile(varargin)

if nargin<1|nargin>29; 
    display('Error in argument list')
end
if nargin>=1
    plotdir = ['../../Manuals/FDS_5_Validation_Guide/FIGURES/'];
    plot_file   = varargin{1};
    devc_col    = varargin{2};
    error_frac  = varargin{3};
    xaxis_label = varargin{4};
    yaxis_label = varargin{5};
    legend_pos  = varargin{6};
    rmin        = varargin{7};
    rmax        = varargin{8};
    dr          = varargin{9};
    xmin        = varargin{10};
    xmax        = varargin{11};
    dx          = varargin{12};
    ymin        = varargin{13};
    ymax        = varargin{14};
    dy          = varargin{15};
    exp_file    = varargin{16};
    exp_format  = varargin{17};
    exp_label   = varargin{18};
    if nargin>18
        fds_file1   = varargin{19};
        fds_format1 = varargin{20};
        fds_label1  = varargin{21};
        nfds        = 1;
    end
    if nargin>21
        fds_file2   = varargin{22};
        fds_format2 = varargin{23};
        fds_label2  = varargin{24};
        nfds        = 2;
    end
    if nargin>24
        fds_file3   = varargin{25};
        fds_format3 = varargin{26};
        fds_label3  = varargin{27};
        nfds        = 3;
    end
    if nargin>27
        fds_file4   = varargin{28};
        fds_format4 = varargin{29};
        fds_label4  = varargin{30};
        nfds        = 4;
    end
end

% experimental data
M = csvread(exp_file,1,0);
r = M(:,1);
w = M(:,2);
e = error_frac*abs(w); % percent error matches DesJardin paper
H(1)=errorbar(r,w,e,exp_format);

% FDS data
if nfds>=1
    hold on
    M = csvread(fds_file1,2,0);
    R = rmin:dr:rmax;
    W = M(3,devc_col);
    H(2)=plot(R,W,fds_format1,'LineWidth',2);
    hold off
end
if nfds>=2
    hold on
    M = csvread(fds_file2,2,0);
    R = rmin:dr:rmax;
    W = M(3,devc_col);
    H(3)=plot(R,W,fds_format2,'LineWidth',2);
    hold off
end
if nfds>=3
    hold on
    M = csvread(fds_file3,2,0);
    R = rmin:dr:rmax;
    W = M(3,devc_col);
    H(4)=plot(R,W,fds_format3,'LineWidth',2);
    hold off
end
if nfds>=4
    hold on
    M = csvread(fds_file4,2,0);
    R = rmin:dr:rmax;
    W = M(3,devc_col);
    H(5)=plot(R,W,fds_format4,'LineWidth',2);
    hold off
end

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

if nfds==1; h = legend(H,exp_label,fds_label1,'Location',legend_pos); end
if nfds==2; h = legend(H,exp_label,fds_label1,fds_label2,'Location',legend_pos); end
if nfds==3; h = legend(H,exp_label,fds_label1,fds_label2,fds_label3,'Location',legend_pos); end
if nfds==4; h = legend(H,exp_label,fds_label1,fds_label2,fds_label3,fds_label4,'Location',legend_pos); end    
set(h,'Interpreter','LaTeX')

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plotdir,plot_file])

