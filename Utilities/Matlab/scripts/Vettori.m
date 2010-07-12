%Myers
%7-06-2010
%Vettori.m

close all
clear all
clc


addpath('../../Validation/Vettori_Flat_Ceiling/Experimental_Data');

%Load FDS Output
M= importdata('Vettori_filenames.csv');
for k = 1:45
fds_data(:,k)=vettori_function(char((M(k))));
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
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

    xmin = 0 ;
    ymin = 0;
    xmax = 220;
    ymax = xmax;
    plot([xmin xmax],[ymin ymax],'k-.');
    axis([xmin xmax ymin ymax]);
    
    xlabel('Measured Activation Time')
    ylabel('Predicted Activation Time')
    
    h =legend(K,'Obstructed Ceiling',' Smooth Ceiling','Location','SouthEast');
    set(h,'Interpreter','LaTeX')
    
    title({'Vettori Flat Ceiling';'Activation Times: Measured v. Predicted'})
    
    hold off
 
    
   
%Print the Plot 
        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[Paper_Width Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);      
        print(gcf,'-dpdf',['../../Manuals/FDS_5_Validation_Guide/FIGURES/Vettori_Flat_Ceiling/Vettori_Activation_Time'])

%Brag about accomplishments
     disp('Great Success')