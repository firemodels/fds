%Myers
%7-06-2010
%vettori_sloped.m

clear all
close all

addpath('../../Validation/Vettori_Sloped_Ceiling');

%Import sorting information
M = importdata('../../Validation/Vettori_Sloped_Ceiling/Experimental_Data/Vettori_Sloped_Guide.csv');

B = M.textdata;

for n=1:72
    %Name FDS results file
    file = ['../../Validation/Vettori_Sloped_Ceiling/FDS_Output_Files/Vettori_Sloped_',char(B(n+1,14)),'.out'];
    
    %Import FDS Data
    %An error will occur at this point if within the simulated time the
    %sprinkler never activated. It is advised that excess time be provided to
    %compensate for this issue and ensure that all sprinklers activate.
    fds_data(:,n) = vettori_sloped_function(file);
    
    %Import Experimental Data
    E(n,:)=str2num(char(B(n+1,8:11)));
end

%Orient Experimental Data
exp_data = E';

%Plot Data, subdivided by obstruction and ceiling slope
for n=1:72
    if n<12
        K(1) = plot(exp_data(1:4,n),fds_data(1:4,n),'r^');
    elseif n>12 & n<=24
        K(2) = plot(exp_data(1:4,n),fds_data(1:4,n),'rd');   
    elseif n>24 & n<=36
        K(3) = plot(exp_data(1:4,n),fds_data(1:4,n),'go');
    elseif n>36 & n<=48
        K(4) = plot(exp_data(1:4,n),fds_data(1:4,n),'gs');
    elseif n>48 & n<60
        K(5) = plot(exp_data(1:4,n),fds_data(1:4,n),'bp');
    elseif n>60 & n<=72
        K(6) = plot(exp_data(1:4,n),fds_data(1:4,n),'bh');
    end
    hold on
end

%Format the Plot 

plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontSize',Title_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

xmin = 0;
ymin = 0;
xmax = 250;
ymax = xmax;
plot([xmin xmax],[ymin ymax],'k-')
axis([xmin xmax ymin ymax])

xlabel('Measured Activation Time (s)','FontSize',Label_Font_Size,'FontName',Font_Name)
ylabel('Predicted Activation Time (s)','FontSize',Label_Font_Size,'FontName',Font_Name)

h = legend(K,'Flat Smooth','Flat Obstructed','13 Smooth','13 Obstructed','24 Smooth','24 Obstructed','Location','SouthEast');

text(0.05*xmax,0.95*ymax,'Vettori Sloped Ceiling Activation Times','FontSize',Title_Font_Size,'FontName',Font_Name)
 
%Print the Plot 

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
display(['Printing plot Vettori_Sloped_Activation_Time ...'])
print(gcf,'-dpdf',['../../Manuals/FDS_Validation_Guide/FIGURES/Vettori_Sloped_Ceiling/Vettori_Sloped_Activation_Time'])

%Brag about accomplishments

disp('Taylor is the Man')
