%Myers
%7-06-2010
%vettori_flat.m

close all
clear all

addpath('../../Validation/Vettori_Flat_Ceiling/Experimental_Data');

%Load FDS Output

M= importdata('Vettori_filenames.csv');
for k = 1:45
    fds_data(:,k)=vettori_flat_function(char((M(k))));
end

for k= 1:45
    %Load Experimental Data
    [exp_header,exp_data] = dvcread('Activation_Times.csv',1);
    %Plot FDS v. Exp. Data
    %Obstructed ceiling cases get blue triangles, smooth get red circles
    
    if k==4||k==5||k==9||k==10||k==14||k==15||k==19||k==20||k==24||k==25||k==29||k==30||k==34||k==35||k==39||k==40||k==44||k==45
        K(1) = plot(exp_data(1:4,k),fds_data(1:4,k),'b^');
    else
        K(2) = plot(exp_data(1:4,k),fds_data(1:4,k),'ro');
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

h = legend(K,'Obstructed Ceiling','Smooth Ceiling','Location','SouthEast');

text(0.05*xmax,0.95*ymax,'Vettori Flat Ceiling Activation Times','FontSize',Title_Font_Size,'FontName',Font_Name)
 
%Print the Plot 

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Width]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Width]);
display(['Printing plot Vettori_Activation_Time ...'])
print(gcf,'-dpdf',['../../Manuals/FDS_Validation_Guide/FIGURES/Vettori_Flat_Ceiling/Vettori_Activation_Time'])

%Brag about accomplishments

disp('Great Success')
