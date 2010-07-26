%Myers
%7-06-2010
%vettori_sloped.m

clear all; close all; clc;
addpath('../../Validation/Vettori_Sloped_Ceiling');

%Import sorting information
[A,B,C] =xlsread('../../Validation/Vettori_Sloped_Ceiling/Experimental_Data/Vettori_Sloped_Guide.csv');

for n=1:72
%Name FDS results file
file = ['../../Validation/Vettori_Sloped_Ceiling/FDS_Output_Files/Vettori_Sloped_',char(B(n+1,14)),'.out'];

%Import FDS Data
%An error will occur at this point if within the simulated time the
%sprinkler never activated. It is advised that excess time be provided to
%compensate for this issue and ensure that all sprinklers activate.
fds_data(:,n) = vettori_sloped_function(file);

%Import Experimental Data
E(n,:)=A(n,8:11);
end

%Orient Experimental Data
exp_data = E';

%Plot Data, subdivided by obstruction and ceiling slope
for n=1:72
K(1) = plot(exp_data(1:4,n),fds_data(1:4,n),'r^');
hold on
end

%Format the Plot 

plot_style
%set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontSize',Title_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[1,1,4.5,4.5])

xmin = 0 ;
ymin = 0;
xmax = 250;
ymax = xmax;
plot([xmin xmax],[ymin ymax],'k-');
axis([xmin xmax ymin ymax]);

xlabel('Measured Activation Time')
ylabel('Predicted Activation Time')

text(0.05*xmax,0.95*ymax,'Vettori Sloped Ceiling Activation Times','FontSize',14,'FontName','Times','Interpreter','LaTeX')

hold off
 
%Print the Plot 

paper_width  = 6.0; % inches
paper_height = 6.0; % inches

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[paper_width paper_height]);
set(gcf,'PaperPosition',[0 0 paper_width paper_height]);      
print(gcf,'-dpdf',['../../Manuals/FDS_5_Validation_Guide/FIGURES/Vettori_Sloped_Ceiling/Vettori_Sloped_Activation_Time'])

%Brag about accomplishments

disp('Vettori Sloped Printed Successfully')