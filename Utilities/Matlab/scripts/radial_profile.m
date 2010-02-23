% McDermott
% 2-15-10
% radial_profile.m
%
% Plot radial profile from experiment vs. FDS results (up to 4)
%
% Example:
%
% function []=radial_profile(plot_file,devc_col,error_frac, ...
%                            xaxis_label,yaxis_label,title_label,text_label,legend_pos, ...
%                            rmin,rmax,dr,xmin,xmax,dx,ymin,ymax,dy, ...
%                            exp_file,exp_format,exp_label, ...
%                            fds_file1,fds_format1,fds_label1, ...
%                            fds_file2,fds_format2,fds_label2, ...
%                            fds_file3,fds_format3,fds_label3, ...
%                            fds_file4,fds_format4,fds_label4)

function []=radial_profile(varargin)

if nargin<1|nargin>32; 
    display('Error in argument list')
end
if nargin>=1
    plotdir = ['../../Manuals/FDS_5_Validation_Guide/FIGURES/'];
    plot_file   = varargin{1};
    devc_col    = varargin{2};
    error_frac  = varargin{3};
    xaxis_label = varargin{4};
    yaxis_label = varargin{5};
    title_label = varargin{6};
    text_label  = varargin{7};
    legend_pos  = varargin{8};
    rmin        = varargin{9};
    rmax        = varargin{10};
    dr          = varargin{11};
    xmin        = varargin{12};
    xmax        = varargin{13};
    dx          = varargin{14};
    ymin        = varargin{15};
    ymax        = varargin{16};
    dy          = varargin{17};
    exp_file    = varargin{18};
    exp_format  = varargin{19};
    exp_label   = varargin{20};
    if nargin>20
        fds_file1   = varargin{21};
        fds_format1 = varargin{22};
        fds_label1  = varargin{23};
        nfds        = 1;
    end
    if nargin>23
        fds_file2   = varargin{24};
        fds_format2 = varargin{25};
        fds_label2  = varargin{26};
        nfds        = 2;
    end
    if nargin>26
        fds_file3   = varargin{27};
        fds_format3 = varargin{28};
        fds_label3  = varargin{29};
        nfds        = 3;
    end
    if nargin>29
        fds_file4   = varargin{30};
        fds_format4 = varargin{31};
        fds_label4  = varargin{32};
        nfds        = 4;
    end
end

% experimental data
M = csvread(exp_file,1,0);
r = M(:,1);
w = M(:,2);
e = error_frac*abs(w); % percent error matches DesJardin paper
e(2:2:length(e)) = 0;  % remove every other error bar for clarity
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

xt = xmin + .05*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,title_label,'FontSize',14,'Interpreter','LaTeX')
xt = xmin + .05*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,text_label,'FontSize',14,'Interpreter','LaTeX')

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel(xaxis_label)
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
legend boxoff

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plotdir,plot_file])

