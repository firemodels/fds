% R. McDermott and C. Cruz
% 6-02-2009
% read_dline.m
%
% Reads verification data configuration file.
%
% Dependencies:
%    define_drow_variables.m
%    verification_data_config_matlab.csv

close all
clear all

vdir = '../../../Verification/';

A = importdata('../verification_data_config_matlab.csv');
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}';

P = textscan(A{2},'%q','delimiter',',');
parameters = P{:}';


if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
    
    define_drow_variables
    
    M = csvread(d1_Filename,2,0);   
    plot(M(:,1),M(:,2))
end