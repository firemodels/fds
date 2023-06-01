% McGrattan 
% 3-12-2021
% tunnel_fire.m
% Test of pressure drop for a tunnel fire
%

close all
clear all

datadir = '../../../out/Moody_Chart/';
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

plot_style

chid='tunnel_fire';
res={'10','20','40'};
lines={'k:','k--','k-'};

for j=1:length(res)
   M = importdata([datadir,chid,'_',res{j},'_line.csv'],',',2);
   x = M.data(:,find(strcmp(M.colheaders,'p-x')));
   p = M.data(:,find(strcmp(M.colheaders,'p')))';
   hh(j)=plot(x,p,lines{j}); hold on;
end

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(hh,'1 m','0.5 m','0.25 m');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

Min_Ind = 0.;
Max_Ind = 1600.;
Min_Dep = 0.;
Max_Dep = 40.;
Title_Position = [0.05 0.90];
X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
text(X_Title_Position,Y_Title_Position,'1600 m Tunnel, 50 MW Fire','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter);

axis([Min_Ind Max_Ind Min_Dep Max_Dep])
xlabel('Distance (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Gauge Pressure (Pa)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(gca,'XTick',[0 400 800 1200 1600])

% add version string if file is available

Git_Filename = [datadir,chid,'_10_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,chid])

