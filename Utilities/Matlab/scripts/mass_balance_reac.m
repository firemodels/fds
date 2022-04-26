% McDermott
% 10 July 2020
% mass_balance_reac.m

close all
clear all

plot_mass_balance('mass_balance_reac','Propane Mass Balance','PROPANE','C3H8');
plot_mass_balance('mass_balance_reac','Oxygen Mass Balance','OXYGEN','O2');
plot_mass_balance('mass_balance_reac','Nitrogen Mass Balance','NITROGEN','N2');
plot_mass_balance('mass_balance_reac','Carbon Dioxide Mass Balance','CARBON DIOXIDE','CO2');
plot_mass_balance('mass_balance_reac','Water Vapor Mass Balance','WATER VAPOR','H2O');
plot_mass_balance('mass_balance_reac','Soot Mass Balance','SOOT','Soot');

function [] = plot_mass_balance(chid,title_text,mass_id,devc_id)

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

ddir = '../../Verification/Species/';
M = importdata([ddir,chid,'_mass.csv'],',',2);

t = M.data(:,1);
m = M.data(:,find(strcmp(M.colheaders,mass_id)));
dmdt = zeros(length(t),1);
for i=2:length(t)
    dmdt(i) = (m(i)-m(i-1))/(t(i)-t(i-1));
end

F = importdata([ddir,chid,'_devc.csv'],',',2);

mdot_out = F.data(:,find(strcmp(F.colheaders,['"',devc_id,' xmax"']))) ...
         + F.data(:,find(strcmp(F.colheaders,['"',devc_id,' xmin"']))) ...
         + F.data(:,find(strcmp(F.colheaders,['"',devc_id,' ymin"']))) ...
         + F.data(:,find(strcmp(F.colheaders,['"',devc_id,' ymax"']))) ...
         + F.data(:,find(strcmp(F.colheaders,['"',devc_id,' Burner"'])));

gen = F.data(:,find(strcmp(F.colheaders,['"',devc_id,' mdot reac"'])));

bal = -dmdt + mdot_out + gen;

%plot(t,zeros(1,length(t)),'k-'); hold on
H(1)=plot(t,dmdt,'g-'); hold on
H(2)=plot(t,mdot_out,'b-');
H(3)=plot(t,gen,'r-');
H(4)=plot(t,bal,'k-');

ylabel('Mass Flow (kg/s)', 'FontSize',Label_Font_Size)
xlabel('Time (s)', 'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'accumulation','in - out','generation', 'balance','location','southwest');
set(lh,'FontSize',Key_Font_Size)

xl=xlim;
yl=ylim;

text(xl(1)+0.05*(xl(2)-xl(1)),yl(1)+0.9*(yl(2)-yl(1)),title_text,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

% check balance and report error

mass_error = abs(mean(bal))/m(end); % error normalized by the total mass in the domain
if mass_error > 1.5e-3
    disp(['Matlab Warning: mass error = ',num2str(mass_error),' in ',chid,' for species ',mass_id])
end

% add version string if file is available

Git_Filename = [ddir,chid,'_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/',chid,'_',devc_id])

end
