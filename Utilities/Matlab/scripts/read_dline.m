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

for i=2:length(A)
    P = textscan(A{i},'%q','delimiter',',');
    parameters = P{:}';
    
    if strcmp(parameters(find(strcmp(headers,'switch_id'))),'d')
        
        define_drow_variables
        
        [H M] = dvcread(d1_Filename);
        d1_Ind_Col = find(strcmp(H,d1_Ind_Col_Name));
        S1 = parse(d1_Dep_Col_Name);
        style = parse(d1_Style);
        for j=1:length(S1)
            d1_Dep_Col = find(strcmp(H,S1(j)));
            if Plot_Type=='linear'
                K(j) = plot(M(:,d1_Ind_Col)/Scale_Ind,M(:,d1_Dep_Col)/Scale_Dep,char(style(j))); hold on
            elseif Plot_Type=='loglog'
                K(j) = loglog(M(:,d1_Ind_Col)/Scale_Ind,M(:,d1_Dep_Col)/Scale_Dep,char(style(j))); hold on
            end
        end
        
        [H M] = dvcread(d2_Filename);
        d2_Ind_Col = find(strcmp(H,d2_Ind_Col_Name));
        S2 = parse(d2_Dep_Col_Name);
        style = parse(d2_Style);
        for j=1:length(S2)
            d2_Dep_Col = find(strcmp(H,S2(j)));
            if Plot_Type=='linear'
                K(length(S1)+j) = plot(M(:,d2_Ind_Col)/Scale_Ind,M(:,d2_Dep_Col)/Scale_Dep,char(style(j)));
            elseif Plot_Type=='loglog'
                K(length(S1)+j) = loglog(M(:,d2_Ind_Col)/Scale_Ind,M(:,d2_Dep_Col)/Scale_Dep,char(style(j)));
            end
        end
        hold off
        
        xlabel(Ind_Title,'interpreter','latex')
        ylabel(Dep_Title,'interpreter','latex')
        axis([Min_Ind Max_Ind Min_Dep Max_Dep])
        if size(Key_Position)>0
            legend(K,[parse(d1_Key),parse(d2_Key)],'Location',Key_Position,'interpreter','latex')
        end
    end
    clear S1 S2 K style H M
    pause
end
