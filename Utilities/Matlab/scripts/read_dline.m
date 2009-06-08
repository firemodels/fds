% R. McDermott and C. Cruz
% 6-02-2009
% read_dline.m
%
% Reads verification data configuration file.
%
% Dependencies:
%    define_drow_variables.m
%    ../verification_data_config_matlab.csv
%    ../functions/dvcread.m
%    ../functions/parse.m

close all
clear all

addpath('../functions')

A = importdata('../verification_data_config_matlab.csv');
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

%for i=2:length(A)
    P = textscan(A{4},'%q','delimiter',',');
    parameters = P{:}';
    
    if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
        
        define_drow_variables
        
        [H M] = dvcread(d1_Filename);
        d1_Ind_Col = find(strcmp(H,d1_Ind_Col_Name));
        S = parse(d1_Dep_Col_Name);
        for j=1:length(S)
            d1_Dep_Col = find(strcmp(H,S(j)));
            K(1) = plot(M(:,d1_Ind_Col)/Scale_Ind,M(:,d1_Dep_Col)/Scale_Dep,'-'); hold on
        end
        
        [H M] = dvcread(d2_Filename);
        d2_Ind_Col = find(strcmp(H,d2_Ind_Col_Name));
        S = parse(d2_Dep_Col_Name);
        for j=1:length(S)
            d2_Dep_Col = find(strcmp(H,S(j)));
            K(2) = plot(M(:,d2_Ind_Col)/Scale_Ind,M(:,d2_Dep_Col)/Scale_Dep,'.:');
        end
        hold off
        
        xlabel(Ind_Title)
        ylabel(Dep_Title)
        axis([Min_Ind Max_Ind Min_Dep Max_Dep])
        legend(K,d1_Key,d2_Key,'Location',Key_Position)
    end
%end