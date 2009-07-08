% McGrattan
% 7-8-2009
% master_validation_script.m
%
% This script processes the plots that are specified in the file 
% validation_data_config_matlab.csv. 

[saved_data,drange] = dataplot('validation',[2:1055]);
scatplot(saved_data,drange,3:40)


