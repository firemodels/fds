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

addpath('../functions')

A = importdata('../verification_data_config_matlab.csv');
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

%for i=2:length(A)
    P = textscan(A{3},'%q','delimiter',',');
    parameters = P{:}';
    
    if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
        
        define_drow_variables
        
        [H M] = dvcread(d1_Filename);
        d1_Ind_Col = find(strcmp(H,d1_Ind_Col_Name));
        d1_Dep_Col = find(strcmp(H,d1_Dep_Col_Name));
        K(1) = plot(M(:,d1_Ind_Col),M(:,d1_Dep_Col),'-'); hold on
        clear H M
        
        [H M] = dvcread(d2_Filename);
        d2_Ind_Col = find(strcmp(H,d2_Ind_Col_Name));
        d2_Dep_Col = find(strcmp(H,d2_Dep_Col_Name));
        K(2) = plot(M(:,d2_Ind_Col),M(:,d2_Dep_Col),'.:');
        hold off
        
        xlabel(Ind_Title)
        ylabel(Dep_Title)
        axis([Min_Ind Max_Ind Min_Dep Max_Dep])
        legend(K,d1_Key,d2_Key,'Location',Key_Position)
    end
%    pause
%end