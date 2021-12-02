% McDermott
% 2 Dec 2021
% mass_balance_gas_volume.m

close all
clear all

plot_mass_balance('mass_balance_gas_volume','');

function [] = plot_mass_balance(chid,title_text)

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

ddir = '../../Verification/Species/';
M = importdata([ddir,chid,'_devc.csv'],',',2);

t = M.data(:,1);
m = M.data(:,find(strcmp(M.colheaders,'"M"')));
dmdt = zeros(length(t),1);
for i=2:length(t)
    dmdt(i) = (m(i)-m(i-1))/(t(i)-t(i-1));
end

mf_x1  = M.data(:,find(strcmp(M.colheaders,'"TMF_X1"')));
mf_x2  = M.data(:,find(strcmp(M.colheaders,'"TMF_X2"')));
mf_y1  = M.data(:,find(strcmp(M.colheaders,'"TMF_Y1"')));
mf_y2  = M.data(:,find(strcmp(M.colheaders,'"TMF_Y2"')));
mf_z1  = M.data(:,find(strcmp(M.colheaders,'"TMF_Z1"')));
mf_z2  = M.data(:,find(strcmp(M.colheaders,'"TMF_Z2"')));

bal = dmdt + (mf_x2-mf_x1 + mf_y2-mf_y1 + mf_z2-mf_z1);

plot(t,zeros(1,length(t)),'k-'); hold on
H(1)=plot(t,dmdt);
H(2)=plot(t,(mf_x2-mf_x1 + mf_y2-mf_y1 + mf_z2-mf_z1));
H(3)=plot(t,bal);

ylabel('mass flow (kg/s)', 'FontSize',Label_Font_Size)
xlabel('time (s)', 'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'dm/dt','out-in','dm/dt+out-in');
set(lh,'FontSize',Key_Font_Size)

text(100,18e-3,title_text,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

% check balance and report error

mass_error = max(abs(bal));
if mass_error > 1e-6
    disp(['Matlab Warning: mass error = ',num2str(mass_error),' in ',chid])
end

% add version string if file is available

Git_Filename = [ddir,chid,'_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',chid])

end
