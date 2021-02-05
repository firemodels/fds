% McDermott
% 12-1-2020
% csiro_flame_position.m
%
% NOTE: this script depends on slread.m, double check "endianness" in slread.m!

close all
clear all

plot_style

expdir = '../../../../exp/CSIRO_Grassland_Fires/';
outdir = '../../../../fds/Validation/CSIRO_Grassland_Fires/Current_Results/';
pltdir = '../../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CSIRO_Grassland_Fires/';

csiro_mesh_extents % returns M1 for Case_F19 and M2 for Case_C064

figure

T = [56 86 138]-0;
for k=1:length(T)
    for n=1:NM3
        nsf = 5;
        % % Case_F19
        % if n>=19 & n<=24
        %     nsf = 10;
        % end
        % Case_F19_fine
        if n>=73 & n<=84
            nsf = 10;
        end
        [Qdata,Time]=slread([outdir 'Case_F19_fine_' num2str(n) '_' num2str(nsf) '.sf'],T(k),T(k)+0.1);
        contourf(M3(n).X,M3(n).Y,Qdata(:,:,1),[20 100]); hold on
    end
end

M = importdata([expdir 'Case_F19_flame_position_data.csv'],',',1);

fl_56_x = M.data(:,find(strcmp(M.colheaders,'56 x')));
fl_56_y = M.data(:,find(strcmp(M.colheaders,'56 y')));
fl_86_x = M.data(:,find(strcmp(M.colheaders,'86 x')));
fl_86_y = M.data(:,find(strcmp(M.colheaders,'86 y')));
fl_138_x = M.data(:,find(strcmp(M.colheaders,'138 x')));
fl_138_y = M.data(:,find(strcmp(M.colheaders,'138 y')));

h(1)=plot(fl_56_x,fl_56_y,'ko','markersize',6); hold on
h(2)=plot(fl_86_x,fl_86_y,'kv','markersize',6);
h(3)=plot(fl_138_x,fl_138_y,'ksq','markersize',8);

plot([0,200],[-100,-100],'k:')
plot([0,200],[+100,+100],'k:')
plot([0,0],[-100,+100],'k:')
plot([200,200],[-100,+100],'k:')

set(gca,'FontName',Font_Name)

set(h(1), 'markerfacecolor', get(h(1), 'color'));
set(h(2), 'markerfacecolor', get(h(2), 'color'));
set(h(3), 'markerfacecolor', get(h(3), 'color'));

axis equal
axis([-10 210 -110 110])
xticks([0:50:200])
yticks([-100:50:100])
xlabel('x (m)','Interpreter',Font_Interpreter)
ylabel('y (m)','Interpreter',Font_Interpreter)

lh=legend(h,'56 s','86 s','138 s','location','southeast');

Paper_Units     = get(gcf,'paperunits');
Paper_Pos       = get(gcf,'paperposition');
Paper_Width     = Paper_Pos(4);
Paper_Height    = Paper_Pos(4);
set(gcf,'Visible','on');
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[pltdir,'Case_F19_flame_position'])

figure

T = [27 53 100];
for k=1:length(T)
    for n=1:NM4
        nsf = 5;
        if n>=19 & n<=24
            nsf = 10;
        end
        [Qdata,Time]=slread([outdir 'Case_C064_fine_' num2str(n) '_' num2str(nsf) '.sf'],T(k),T(k)+0.1);
        contourf(M4(n).X,M4(n).Y,Qdata(:,:,1),[20 100]); hold on
    end
end

M = importdata([expdir 'Case_C064_flame_position_data.csv'],',',1);

fl_27_x = M.data(:,find(strcmp(M.colheaders,'27 x')));
fl_27_y = M.data(:,find(strcmp(M.colheaders,'27 y')));
fl_53_x = M.data(:,find(strcmp(M.colheaders,'53 x')));
fl_53_y = M.data(:,find(strcmp(M.colheaders,'53 y')));
fl_100_x = M.data(:,find(strcmp(M.colheaders,'100 x')));
fl_100_y = M.data(:,find(strcmp(M.colheaders,'100 y')));

h(1)=plot(fl_27_x,fl_27_y,'ko','markersize',6); hold on
h(2)=plot(fl_53_x,fl_53_y,'kv','markersize',6);
h(3)=plot(fl_100_x,fl_100_y,'ksq','markersize',8);

plot([0,100],[-50,-50],'k:')
plot([0,100],[+50,+50],'k:')
plot([0,0],[-50,+50],'k:')
plot([100,100],[-50,+50],'k:')

set(gca,'FontName',Font_Name)

set(h(1), 'markerfacecolor', get(h(1), 'color'));
set(h(2), 'markerfacecolor', get(h(2), 'color'));
set(h(3), 'markerfacecolor', get(h(3), 'color'));

axis equal
axis([-5 105 -55 55])
xticks([0:10:100])
yticks([-50:10:50])
xlabel('x (m)','Interpreter',Font_Interpreter)
ylabel('y (m)','Interpreter',Font_Interpreter)

lh=legend(h,'27 s','53 s','100 s','location','southeast');

Paper_Units     = get(gcf,'paperunits');
Paper_Pos       = get(gcf,'paperposition');
Paper_Width     = Paper_Pos(4);
Paper_Height    = Paper_Pos(4);
set(gcf,'Visible','on');
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

print(gcf,'-dpdf',[pltdir,'Case_C064_flame_position'])





