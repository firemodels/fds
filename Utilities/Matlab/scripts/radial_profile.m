% McDermott
% 2-15-10
% radial_profile.m
%
% Plot radial profile from experiment vs. FDS results (up to 5)
%
% Example:
%
% function []=radial_profile(plot_file,data_format,devc_col,error,error_type,tmin, ...
%                            xaxis_label,yaxis_label,title_label,text_label,legend_pos, ...
%                            rmin,rmax,dr,xmin,xmax,dx,ymin,ymax,dy,git_file, ...
%                            exp_file,exp_format,exp_label,exp_stride, ...
%                            fds_file1,fds_format1,fds_label1, ...
%                            fds_file2,fds_format2,fds_label2, ...
%                            fds_file3,fds_format3,fds_label3, ...
%                            fds_file4,fds_format4,fds_label4, ...
%                            fds_file5,fds_format5,fds_label5)

function []=radial_profile(varargin)

if nargin<1|nargin>100;
    display('Error in argument list')
end
if nargin>=1
    iarg = 1;
    plot_file   = varargin{iarg}; iarg=iarg+1;
    data_format = varargin{iarg}; iarg=iarg+1;
    devc_col    = varargin{iarg}; iarg=iarg+1;
    error       = varargin{iarg}; iarg=iarg+1;
    error_type  = varargin{iarg}; iarg=iarg+1;
    tmin        = varargin{iarg}; iarg=iarg+1;
    xaxis_label = varargin{iarg}; iarg=iarg+1;
    yaxis_label = varargin{iarg}; iarg=iarg+1;
    title_label = varargin{iarg}; iarg=iarg+1;
    text_label  = varargin{iarg}; iarg=iarg+1;
    legend_pos  = varargin{iarg}; iarg=iarg+1;
    rmin        = varargin{iarg}; iarg=iarg+1;
    rmax        = varargin{iarg}; iarg=iarg+1;
    dr          = varargin{iarg}; iarg=iarg+1;
    xmin        = varargin{iarg}; iarg=iarg+1;
    xmax        = varargin{iarg}; iarg=iarg+1;
    dx          = varargin{iarg}; iarg=iarg+1;
    ymin        = varargin{iarg}; iarg=iarg+1;
    ymax        = varargin{iarg}; iarg=iarg+1;
    dy          = varargin{iarg}; iarg=iarg+1;
    git_file    = varargin{iarg}; iarg=iarg+1;
    exp_file    = varargin{iarg}; iarg=iarg+1;
    exp_format  = varargin{iarg}; iarg=iarg+1;
    exp_label   = varargin{iarg}; iarg=iarg+1;
    exp_stride  = varargin{iarg}; iarg=iarg+1;
    nfds        = 0;
    if nargin>iarg
        fds_file1   = varargin{iarg}; iarg=iarg+1;
        fds_format1 = varargin{iarg}; iarg=iarg+1;
        fds_label1  = varargin{iarg}; iarg=iarg+1;
        nfds        = 1;
    end
    if nargin>iarg
        fds_file2   = varargin{iarg}; iarg=iarg+1;
        fds_format2 = varargin{iarg}; iarg=iarg+1;
        fds_label2  = varargin{iarg}; iarg=iarg+1;
        nfds        = 2;
    end
    if nargin>iarg
        fds_file3   = varargin{iarg}; iarg=iarg+1;
        fds_format3 = varargin{iarg}; iarg=iarg+1;
        fds_label3  = varargin{iarg}; iarg=iarg+1;
        nfds        = 3;
    end
    if nargin>iarg
        fds_file4   = varargin{iarg}; iarg=iarg+1;
        fds_format4 = varargin{iarg}; iarg=iarg+1;
        fds_label4  = varargin{iarg}; iarg=iarg+1;
        nfds        = 4;
    end
    if nargin>iarg
        fds_file5   = varargin{iarg}; iarg=iarg+1;
        fds_format5 = varargin{iarg}; iarg=iarg+1;
        fds_label5  = varargin{iarg}; iarg=iarg+1;
        nfds        = 5;
    end
end

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% experimental data
M = csvread(exp_file,1,0);
r = M(1:exp_stride:end,1);
w = M(1:exp_stride:end,2);
if (error>0)
    if strcmp(error_type,'relative')
        e = error*abs(w); % percent error matches DesJardin paper
    else
        e = error*ones(1,length(w)/length(error));
    end
    % put error bars on every other data point for readability
    plot(r,w,exp_format); hold on
    H(1)=errorbar(r,w,-e,+e,exp_format);
else
    H(1)=plot(r,w,exp_format);
end

% FDS data
if nfds>=1
    hold on
    M = csvread(fds_file1,2,0);
    if strcmp(data_format,'row')
        R = rmin:dr:rmax;
        W = M(3,devc_col);
    elseif strcmp(data_format,'col')
        R = M(:,1);
        W = M(:,devc_col);
    elseif strcmp(data_format,'line')
        R = M(:,devc_col(1));
        W = M(:,devc_col(2));
    elseif strcmp(data_format,'mean')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = mean(M(T,devc_col),1);
    elseif strcmp(data_format,'rms')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = std(M(T,devc_col),1);
    elseif strcmp(data_format,'tke')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        urange = devc_col(1):devc_col(2);
        vrange = devc_col(3):devc_col(4);
        u = M(T,urange)-ones(length(T),1)*mean(M(T,urange),1);
        v = M(T,vrange)-ones(length(T),1)*mean(M(T,vrange),1);
        W = mean(0.5*(u.*u+v.*v),1);
    elseif strcmp(data_format,'col_tke')
        R = M(:,1);
        u = M(:,devc_col(1));
        v = M(:,devc_col(2));
        W = 0.5*(u.*u + v.*v);
    end
    H(2)=plot(R,W,fds_format1,'LineWidth',1);
    hold off
end
if nfds>=2
    hold on
    M = csvread(fds_file2,2,0);
    if strcmp(data_format,'row')
        R = rmin:dr:rmax;
        W = M(3,devc_col);
    elseif strcmp(data_format,'col')
        R = M(:,1);
        W = M(:,devc_col);
    elseif strcmp(data_format,'line')
        R = M(:,devc_col(1));
        W = M(:,devc_col(2));
    elseif strcmp(data_format,'mean')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = mean(M(T,devc_col),1);
    elseif strcmp(data_format,'rms')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = std(M(T,devc_col),1);
    elseif strcmp(data_format,'tke')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        urange = devc_col(1):devc_col(2);
        vrange = devc_col(3):devc_col(4);
        u = M(T,urange)-ones(length(T),1)*mean(M(T,urange),1);
        v = M(T,vrange)-ones(length(T),1)*mean(M(T,vrange),1);
        W = mean(0.5*(u.*u+v.*v),1);
    elseif strcmp(data_format,'col_tke')
        R = M(:,1);
        u = M(:,devc_col(1));
        v = M(:,devc_col(2));
        W = 0.5*(u.*u + v.*v);
    end
    H(3)=plot(R,W,fds_format2,'LineWidth',1);
    hold off
end
if nfds>=3
    hold on
    M = csvread(fds_file3,2,0);
    if strcmp(data_format,'row')
        R = rmin:dr:rmax;
        W = M(3,devc_col);
    elseif strcmp(data_format,'col')
        R = M(:,1);
        W = M(:,devc_col);
    elseif strcmp(data_format,'line')
        R = M(:,devc_col(1));
        W = M(:,devc_col(2));
    elseif strcmp(data_format,'mean')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = mean(M(T,devc_col),1);
    elseif strcmp(data_format,'rms')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = std(M(T,devc_col),1);
    elseif strcmp(data_format,'tke')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        urange = devc_col(1):devc_col(2);
        vrange = devc_col(3):devc_col(4);
        u = M(T,urange)-ones(length(T),1)*mean(M(T,urange),1);
        v = M(T,vrange)-ones(length(T),1)*mean(M(T,vrange),1);
        W = mean(0.5*(u.*u+v.*v),1);
    elseif strcmp(data_format,'col_tke')
        R = M(:,1);
        u = M(:,devc_col(1));
        v = M(:,devc_col(2));
        W = 0.5*(u.*u + v.*v);
    end
    H(4)=plot(R,W,fds_format3,'LineWidth',1);
    hold off
end
if nfds>=4
    hold on
    M = csvread(fds_file4,2,0);
    if strcmp(data_format,'row')
        R = rmin:dr:rmax;
        W = M(3,devc_col);
    elseif strcmp(data_format,'col')
        R = M(:,1);
        W = M(:,devc_col);
    elseif strcmp(data_format,'line')
        R = M(:,devc_col(1));
        W = M(:,devc_col(2));
    elseif strcmp(data_format,'mean')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = mean(M(T,devc_col),1);
    elseif strcmp(data_format,'rms')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = std(M(T,devc_col),1);
    elseif strcmp(data_format,'tke')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        urange = devc_col(1):devc_col(2);
        vrange = devc_col(3):devc_col(4);
        u = M(T,urange)-ones(length(T),1)*mean(M(T,urange),1);
        v = M(T,vrange)-ones(length(T),1)*mean(M(T,vrange),1);
        W = mean(0.5*(u.*u+v.*v),1);
    elseif strcmp(data_format,'col_tke')
        R = M(:,1);
        u = M(:,devc_col(1));
        v = M(:,devc_col(2));
        W = 0.5*(u.*u + v.*v);
    end
    H(5)=plot(R,W,fds_format4,'LineWidth',1);
    hold off
end
if nfds>=5
    hold on
    M = csvread(fds_file5,2,0);
    if strcmp(data_format,'row')
        R = rmin:dr:rmax;
        W = M(3,devc_col);
    elseif strcmp(data_format,'col')
        R = M(:,1);
        W = M(:,devc_col);
    elseif strcmp(data_format,'line')
        R = M(:,devc_col(1));
        W = M(:,devc_col(2));
    elseif strcmp(data_format,'mean')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = mean(M(T,devc_col),1);
    elseif strcmp(data_format,'rms')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        W = std(M(T,devc_col),1);
    elseif strcmp(data_format,'tke')
        R = rmin:dr:rmax;
        T = find(M(:,1)>=tmin);
        urange = devc_col(1):devc_col(2);
        vrange = devc_col(3):devc_col(4);
        u = M(T,urange)-ones(length(T),1)*mean(M(T,urange),1);
        v = M(T,vrange)-ones(length(T),1)*mean(M(T,vrange),1);
        W = mean(0.5*(u.*u+v.*v),1);
    elseif strcmp(data_format,'col_tke')
        R = M(:,1);
        u = M(:,devc_col(1));
        v = M(:,devc_col(2));
        W = 0.5*(u.*u + v.*v);
    end
    H(6)=plot(R,W,fds_format5,'LineWidth',1);
    hold off
end

xt = xmin + .05*(xmax-xmin);
yt = ymin + .92*(ymax-ymin);
text(xt,yt,title_label,'FontName',Font_Name,'FontSize',14,'Interpreter',Font_Interpreter)
xt = xmin + .05*(xmax-xmin);
yt = ymin + .84*(ymax-ymin);
text(xt,yt,text_label,'FontName',Font_Name,'FontSize',14,'Interpreter',Font_Interpreter)

xlabel(xaxis_label,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
ylabel(yaxis_label,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
axis([xmin xmax ymin ymax])
set(gca,'YTick',ymin:dy:ymax)
set(gca,'YMinorTick','on')
set(gca,'XTick',xmin:dx:xmax)
set(gca,'XMinorTick','on')
%set(gca,'YMinorGrid','on')

if nfds==0; h = legend(H,exp_label,'Location',legend_pos); end
if nfds==1; h = legend(H,exp_label,fds_label1,'Location',legend_pos); end
if nfds==2; h = legend(H,exp_label,fds_label1,fds_label2,'Location',legend_pos); end
if nfds==3; h = legend(H,exp_label,fds_label1,fds_label2,fds_label3,'Location',legend_pos); end
if nfds==4; h = legend(H,exp_label,fds_label1,fds_label2,fds_label3,fds_label4,'Location',legend_pos); end
if nfds==5; h = legend(H,exp_label,fds_label1,fds_label2,fds_label3,fds_label4,fds_label5,'Location',legend_pos); end
set(h,'Interpreter',Font_Interpreter)
set(h,'FontSize',Key_Font_Size)

% add Git revision if file is available

addverstr(gca,git_file,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',plot_file)

