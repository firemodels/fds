% McGrattan
% 7-8-2009
% master_validation_script.m
%
% This script processes the plots that are specified in the file 
% validation_data_config_matlab.csv.
%
% See also: dataplot, scatplot

close all
clear all

addpath 'scripts'

[saved_data,drange] = dataplot('validation');
scatplot(saved_data,drange)


