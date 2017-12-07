% McDermott
% 05 Dec 2017
% mass_balance.m

close all
clear all

plot_mass_balance('mass_flux_wall_yindex','Primitive Species Mass Balance');
plot_mass_balance('mass_flux_wall_zindex','Lumped Species Mass Balance');


function [] = plot_mass_balance(chid,title_text)

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

ddir = '../../Verification/Flowfields/';
M = importdata([ddir,chid,'_mass.csv'],',',2);

t = M.data(:,1);
m = M.data(:,find(strcmp(M.colheaders,'WATER VAPOR')));
dmdt = zeros(length(t),1);
for i=2:length(t)
    dmdt(i) = (m(i)-m(i-1))/(t(i)-t(i-1));
end

F = importdata([ddir,chid,'_devc.csv'],',',2);

mdot_in  = F.data(:,find(strcmp(F.colheaders,'"H2O in"')));
mdot_out = F.data(:,find(strcmp(F.colheaders,'"H2O out"')));

bal = dmdt - mdot_in - mdot_out;

plot(t,zeros(1,length(t)),'k-'); hold on
H(1)=plot(t,mdot_in);
H(2)=plot(t,-mdot_out);
H(3)=plot(t,bal);

ylabel('mass flow (kg/s)', 'FontSize',Label_Font_Size)
xlabel('time (s)', 'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Inlet H2O','Outlet H2O','dm/dt+out-in');
set(lh,'FontSize',Key_Font_Size)

text(100,18e-3,title_text,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

% check balance and report error

mass_error = abs(mean(bal(find(t>1000))));
if mass_error > 1e-5
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
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',chid,'_mass_balance'])

end
