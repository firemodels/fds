% McDermott
% 6-02-2009
% read_dline.m
%
% Reads verification data configuration file.

close all
clear all

addpath('../functions')
vdir = '../../../Verification/';

A = importdata('../verification_data_config_matlab.csv');
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}';

P = textscan(A{2},'%q','delimiter',',');
parameters = P{:}';


if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
    
    Quantity        = char(parameters(find(strcmp(headers,'Quantity'))));
    Group           = char(parameters(find(strcmp(headers,'Group'))));
    d1_Filename     = [vdir,char(parameters(find(strcmp(headers,'d1_Filename'))))];
    d1_Col_Name_Row = str2num(char(parameters(find(strcmp(headers,'d1_Col_Name_Row')))))
    
    M = csvread(d1_Filename,2,0);
    
    plot(M(:,1),M(:,2))
end