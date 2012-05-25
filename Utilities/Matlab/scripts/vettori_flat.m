%Myers
%7-06-2010
%vettori_flat.m

close all
clear all
clc

addpath('../../Validation/Vettori_Flat_Ceiling/Experimental_Data');

%Load FDS Output

M= importdata('Vettori_filenames.csv');
for k = 1:45
    fds_data(:,k)=vettori_flat_function(char((M(k))));
end

for k= 1:45
    hold on
    %Load Experimental Data
    [exp_header,exp_data] = dvcread('Activation_Times.csv',1);
    %Plot FDS v. Exp. Data
    %Obstructed ceiling cases get blue triangles, smooth get red circles
    
    if k==4||k==5||k==9||k==10||k==14||k==15||k==19||k==20||k==24||k==25||k==29||k==30||k==34||k==35||k==39||k==40||k==44||k==45
        K(1) = plot(exp_data(1:4,k),fds_data(1:4,k),'b^');
    else
        K(2)= plot(exp_data(1:4,k),fds_data(1:4,k),'ro');
    end
end

%Format the Plot 

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontSize',Title_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[1,1,4.5,4.5])

axis on
axis equal

xmin = 0;
ymin = 0;
xmax = 250;
ymax = xmax;
plot([xmin xmax],[ymin ymax],'k-');
axis([xmin xmax ymin ymax]);

xlabel('Measured Activation Time (s)')
ylabel('Predicted Activation Time (s)')

h = legend(K,'Obstructed Ceiling','Smooth Ceiling','Location','SouthEast');
set(h,'Interpreter',Font_Interpreter)

text(0.05*xmax,0.95*ymax,'Vettori Flat Ceiling Activation Times','FontSize',14,'FontName','Times','Interpreter',Font_Interpreter)

hold off
 
%Print the Plot 

paper_width  = 6.0; % inches
paper_height = 6.0; % inches

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[paper_width paper_height]);
set(gcf,'PaperPosition',[0 0 paper_width paper_height]);      
print(gcf,'-dpdf',['../../Manuals/FDS_Validation_Guide/FIGURES/Vettori_Flat_Ceiling/Vettori_Activation_Time'])

%Brag about accomplishments

disp('Great Success')