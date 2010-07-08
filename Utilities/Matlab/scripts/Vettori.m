%Myers
%7-06-2010
%Vettori.m

close all
clear all
clc


addpath('../../../Validation/Vettori_Flat_Ceiling/Experimental_Data');

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
 
    plot(exp_data(1:4,k),fds_data(1:4,k),'rp');
 
end

%Format plot    
    xmin = 0 ;
    ymin = 0;
    xmax = 220;
    ymax = xmax;
    plot([xmin xmax],[ymin ymax],'k-.');
    axis([xmin xmax ymin ymax]);
    
    xlabel('Measured Activation Time')
    ylabel('Predicted Activation Time')
    
    TITLE(['Activation Times: Measured v. Predicted'])
 hold off
 
%Print the Plot 
     print(gcf,'-dpdf',['../../../Manuals/FDS_5_Validation_Guide/FIGURES/Vettori_Flat_Ceiling/Vettori_Activation_Time'])

 %Brag about accomplishments
     disp('Great Success')