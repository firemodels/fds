% $Id$
% $Revision$

% Generate fds summary data for scatterplots

outdir = '../../../out/FM_FPRF_Datacenter/FDS_Output_Files/';
expdir = '../../../exp/FM_FPRF_Datacenter/';

[exp_data] = csvread([expdir,'/fm_exp.csv'],1);

% Low flow test
[fds_data] = csvread([outdir,'FM_Datacenter_Veltest_Low_devc.csv'],13);
n_fds_data=size(fds_data,1);
%compute average pressures
for i=2:7
   fds_avg(i-1) = mean(fds_data(1:n_fds_data,i));
end

fds_header={'Time','Low SF-CA','Low HA-CP','High SF-CA','High HA-CP'};

fds_out(1)=100;
fds_out(2)=fds_avg(1)-fds_avg(2);
fds_out(3)=fds_avg(5)-fds_avg(6);

% High flow test
[fds_data] = csvread([outdir,'FM_Datacenter_Veltest_High_devc.csv'],13);

%compute average pressures
for i=2:7
   fds_avg(i-1) = mean(fds_data(1:n_fds_data,i));
end

fds_out(4)=fds_avg(1)-fds_avg(2);
fds_out(5)=fds_avg(5)-fds_avg(6);

fid=fopen([outdir,'FM_Datacenter_fds_data.csv'],'w+');

fprintf(fid,'%s, %s, %s, %s, %s \n',fds_header{1:5});
fprintf(fid,'%f, %f, %f, %f, %f, \n',fds_out(1:5));

fclose(fid);

%get soot values

[fds_data] = csvread([outdir,'FM_Datacenter_Low_C3H6_SF_devc.csv'],15);

fds_out(6) = mean(fds_data(:,28))*1000000;
fds_out(7) = mean(fds_data(:,46))*1000000;
fds_out(8) = mean(fds_data(:,65))*1000000;

[fds_data] = csvread([outdir,'FM_Datacenter_High_C3H6_SF_devc.csv'],15);

fds_out(9)  = mean(fds_data(:,28))*1000000;
fds_out(10) = mean(fds_data(:,46))*1000000;
fds_out(11) = mean(fds_data(:,65))*1000000;

[fds_data] = csvread([outdir,'FM_Datacenter_Low_C3H6_HA_devc.csv'],15);

fds_out(12) = mean(fds_data(:,28))*1000000;
fds_out(13) = mean(fds_data(:,46))*1000000;
fds_out(14) = mean(fds_data(:,65))*1000000;

[fds_data] = csvread([outdir,'FM_Datacenter_High_C3H6_HA_devc.csv'],15);

fds_out(15) = mean(fds_data(:,28))*1000000;
fds_out(16) = mean(fds_data(:,46))*1000000;
fds_out(17) = mean(fds_data(:,65))*1000000;

[fds_data] = csvread([outdir,'FM_Datacenter_Low_Cable_SF_devc.csv'],13);

fds_out(18) = mean(fds_data(:,28))*1000000;
fds_out(19) = mean(fds_data(:,46))*1000000;
fds_out(20) = mean(fds_data(:,65))*1000000;

[fds_data] = csvread([outdir,'FM_Datacenter_High_Cable_SF_devc.csv'],13);

fds_out(21) = mean(fds_data(:,28))*1000000;
fds_out(22) = mean(fds_data(:,46))*1000000;
fds_out(23) = mean(fds_data(:,65))*1000000;

x=[0.01 0.122 0.2 0.3 0.5 1 5 10 50 100 500 1000];
logx=min(max(-2.7,log(x+0.00001)),1);
errx=0.01*(3.8184*logx.^2-7.7783*logx+14.346);
toterr=(errx.^2+0.1^2+0.1^2+0.05^2).^0.5;
xerrp = x + 2*toterr.*x;
xerrm = max(0.00001,x - 2*toterr.*x);

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

hx=loglog(x,x,'k-',x,xerrp,'k--',x,xerrm,'k--');
hold on;
h1=loglog(exp_data(6:8),fds_out(6:8),'ro');
h2=loglog(exp_data(9:11),fds_out(9:11),'r+');
h3=loglog(exp_data(12:14),fds_out(12:14),'bo');
h4=loglog(exp_data(15:17),fds_out(15:17),'b+');
h5=loglog(exp_data(18:20),fds_out(18:20),'go');
h6=loglog(exp_data(21:23),fds_out(21:23),'g+');

xlim([0.01 300]);
ylim([0.01 300]);

lh=legend([h1 h2 h3 h4 h5 h6],...
    'C3H6 Low SF','C3H6 High SF','C3H6 Low HA','C3H6 High HA','Cable Low SF','Cable High SF',...
    'Location','southeast');
set(lh,'FontSize',Key_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)
xtitle = ['Measured Soot Concentration (mg/m^3)'];
ytitle = ['Predicted Soot Concentration (mg/m^3)'];
xlabel(xtitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel(ytitle,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

git_file = [outdir,'FM_Datacenter_Low_C3H6_SF_git.txt'];
addverstr(gca,git_file,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Soot'];
print(gcf,'-dpdf',plotname);
hold off

filename = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/pressure.tex';
fid = fopen(filename,'wt');

pres_dump(1)=fds_out(2);
pres_dump(2)=exp_data(2);
pres_dump(3)=exp_data(2)*.19;
pres_dump(4)=fds_out(3);
pres_dump(5)=exp_data(3);
pres_dump(6)=exp_data(3)*.19;
pres_dump(7)=fds_out(4);
pres_dump(8)=exp_data(4);
pres_dump(9)=exp_data(4)*.1;
pres_dump(10)=fds_out(5);
pres_dump(11)=exp_data(5);
pres_dump(12)=exp_data(5)*.1;


fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|c|} \hline');
fprintf(fid,'%s\n','Fan Speed & FDS SF to CA & Exp SF to CA & FDS HA to CP & Exp HA to CP \\');
fprintf(fid,'%s\n','  & (Pa) & (Pa) & (Pa) & (Pa)  \\ \hline\hline');
fprintf(fid,'78 ACH & %5.1f & %5.1f $\\pm$ %5.1f & %5.1f & %5.1f $\\pm$ %5.1f %s\n',pres_dump(1:6),'\\');
fprintf(fid,'265 ACH & %5.1f & %5.1f $\\pm$ %5.1f & %5.1f & %5.1f $\\pm$ %5.1f %s\n',pres_dump(7:12),'\\');
fprintf(fid,'%s\n','\hline');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');
