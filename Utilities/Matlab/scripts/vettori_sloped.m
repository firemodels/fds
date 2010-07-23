%Myers
%7-06-2010
%vettori_sloped.m

clear all; close all; clc;
addpath('../../Validation/Vettori_Sloped_Ceiling');

%Import sorting information
[A,B,C] =xlsread('Vettori_Sloped_Guide.csv');


for n=1:72
%Name FDS results file
file = ['../../Validation/Vettori_Sloped_Ceiling/FDS_Output/Vettori_Sloped_',char(B(n+1,14)),'.out'];

%Import FDS Data
fds_data(:,n) = vettori_sloped_function(file);

%Import Experimental Data
E(n,:)=A(n,8:11);
end

%Orient Experimental Data
exp_data = D';

%Plot Data, subdivided by obstruction and ceiling slope
for n=1:72
if n<=12
    K(1) = plot(exp_data(1:4,n),fds_data(1:4,n),'r^');
elseif 12<n<=24
    K(2) = plot(exp_data(1:4,n),fds_data(1:4,n),'rd')
elseif 24<n<=36
    K(3) = plot(exp_data(1:4,n),fds_data(1:4,n),'bo')
elseif 36<n<=48
    K(4) = plot(exp_data(1:4,n),fds_data(1:4,n),'bs')
elseif 48<n<=60
    K(5) = plot(exp_data(1:4,n),fds_data(1:4,n),'gp')
else
    K(3) = plot(exp_data(1:4,n),fds_data(1:4,n),'gh')
end

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

h =legend(K,'Flat Smooth','Flat Obstructed','13 Degree Smooth','13 Degree Obstructed','24 Degree Smooth','24 Degree Obstructed','Location','SouthEast');
set(h,'Interpreter','LaTeX')

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