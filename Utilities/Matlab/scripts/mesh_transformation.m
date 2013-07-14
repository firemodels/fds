% McDermott and McGrattan
% 10-4-12
% mesh_transformation.m
%
% Postprocessing script to extract transformation mapping from SMV.

close all
clear all

% Read SMV file and store transformations

datadir = '../../Verification/Miscellaneous/';
fid  = fopen([datadir,'mesh_transformation.smv'],'r');

I = 50;
J = 50;
LX = 1.5;
LY = 1.5;
dx = LX/I;
dy = LY/J;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if strcmp(tline,'TRNX')
        %disp(tline)
        for i=1:3
            skip=fgetl(fid);
        end
        for i=1:I+1
            tline=fgetl(fid);
            C=textscan(tline,'%f');
            %disp(C{1})
            cx(i) = C{1}(1)*dx;
            x(i) = C{1}(2);
        end
    end
    if strcmp(tline,'TRNY')
        %disp(tline)
        for j=1:3
            skip=fgetl(fid);
        end
        for j=1:J+1
            tline=fgetl(fid);
            C=textscan(tline,'%f');
            %disp(C{1})
            cy(j) = C{1}(1)*dy;
            y(j) = C{1}(2);
        end
    end
end
fclose(fid);

% Make the plots

% From the mesh_transformation.fds input file:
CC = [0.3 1.2];
PC = [0.5 1.0];

plotdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';
plot_style
Title_Font_Size=20;

XMinorTick = 0.1:0.1:1.4;

plot(cx,x,'k-','LineWidth',2); hold on
for i=1:length(XMinorTick)
    % piece-wise linear function
    if XMinorTick(i)<CC(1)
        YMinorTick = (PC(1)/CC(1))*XMinorTick(i);
    elseif XMinorTick(i)<CC(2)
        YMinorTick = PC(1) + (PC(2)-PC(1))/(CC(2)-CC(1))*(XMinorTick(i)-CC(1));
    elseif XMinorTick(i)<=LX
        YMinorTick = PC(2) + (LX-PC(2))/(LX-CC(2))*(XMinorTick(i)-CC(2));
    end
    plot([XMinorTick(i) XMinorTick(i)],[0 YMinorTick],'k-')
    plot([0 XMinorTick(i)],[YMinorTick YMinorTick],'k-')
end
set(gca,'XTick',[0:0.3:1.5])
set(gca,'YTick',[0:0.3:1.5])
xlabel('\xi','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
ylabel('\itx','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

set(gca,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf',[plotdir,'piece_wise_linear_trnx'])

figure

plot(cy,y,'k-','LineWidth',2); hold on
for i=1:length(XMinorTick)
    % polynomial function
    YMinorTick = 2*XMinorTick(i) - 2*XMinorTick(i)^2 + 0.8889*XMinorTick(i)^3;
    plot([XMinorTick(i) XMinorTick(i)],[0 YMinorTick],'k-')
    plot([0 XMinorTick(i)],[YMinorTick YMinorTick],'k-')
end
set(gca,'XTick',[0:0.3:1.5])
set(gca,'YTick',[0:0.3:1.5])
xlabel('\xi','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
ylabel('\itx','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

set(gca,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf',[plotdir,'polynomial_trnx'])




