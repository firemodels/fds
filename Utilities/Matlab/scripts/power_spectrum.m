% McDermott
% 2-26-10
% power_spectrum.m

function []=power_spectrum(spectrum_file,dt_puff,slope_scale,xaxis_title,yaxis_title,xmin,xmax,ymin,ymax, ...
                           devc_file,devc_col,tmin,tmax,dt,spectrum_style,nyquist_style, ...
                           title1,title2,git_file)

M = csvread(devc_file,2,0);

range = find(M(:,1)>tmin & M(:,1)<tmax);

% power spectrum
for i=devc_col
    W = M(range,i);
    n = length(W);
    u = fft(W);
    p = u.*conj(u)/n^2;
    if i==devc_col(1)
        pave=p;
    else
        pave = pave + p;
    end
end
pave = pave/length(devc_col);

t = M(range,1);
T = t(length(t))-t(1);
n2 = floor(n/2)+1;
k = 2:n2;
f = k/T;
f_puff = 1/dt_puff;
dt_out = 19/1000;

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

kk = find(f>5 & f<90);
H(1)=loglog(f,pave(k),spectrum_style); hold on
k_fds = find(pave(k)==max(pave(k)));
f(k_fds);

loglog(f(kk),slope_scale*f(kk).^(-5/3),'k-','LineWidth',1)
loglog([f_puff,f_puff],[min(pave(k)),max(pave(k))],'k--','LineWidth',2);
%loglog([.5/dt,.5/dt],[min(pave(k)),1e-3*max(pave(k))],nyquist_style,'LineWidth',2);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel(xaxis_title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(yaxis_title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([xmin xmax ymin ymax])
%h=legend(H,'FDS W, 1.5 cm','Location','Northeast'); set(h,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
%legend boxoff

xt = 1.5e-1;
yt = 3.5e0;
text(xt,yt,title1,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
xt = 1.5e-1;
yt = .9e0;
text(xt,yt,title2,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)

text(1.5e-1,3.5e-4,'measured','FontSize',Label_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);
text(1.5e-1,1.2e-4,'puffing','FontSize',Label_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);
text(1.5e-1,.4e-4,'frequency','FontSize',Label_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);
annotation('arrow',[.3 .4],[.3 .3]);
text(1e1,4e-2,'-5/3','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
%text(1e2,1.1e-3,'Nyquist','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)
%text(1.4e2,.4e-3,'limit','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision
addverstr(gca,git_file,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',spectrum_file)




