% R. McDermott
% 6-11-2009
% scatterplot.m
%
% Generate scatter plots.  Must first run read_dline.m to generate
% measured_metric and predicted_metric.
%
% Dependencies:
%    define_qrow_variables.m
%    ../scatterplot_config_matlab.csv

close all

qfil = ['../scatterplot_config_matlab.csv'];
qrange = 2:2;

addpath('../functions')

Q = importdata(qfil);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

for j=qrange
    if j>length(Q); break; end
    
    define_qrow_variables
    
    for i=drange
        if strcmp(Save_Quantity(i),Quantity_Label)
            plot(Save_Measured_Metric(i),Save_Predicted_Metric(i),'b^'); hold on
        end
    end
    
end