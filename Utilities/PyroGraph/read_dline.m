% McDermott
% 5-22-2009
% read_dline.m

close all
clear all

[M T] = xlsread('verification_data_config_matlab.csv'); % M = num array, T = cell array

for k = 2:length(T(:,5))
    if strcmp((T(k,1)),'d')
        [D H] = xlsread(['../../Verification/',char(T(k,5))]);
        C = textscan(char(T(k,8)),'%s','delimiter',',');
        C = C{:}';
        
        for i=1:length(C)
            for j=2:length(H)
                if strcmp(C(i),H(j))
                    plot(D(:,1),D(:,j)); hold on
                end
            end
        end

        xlabel(char(T(k,7)))
        ylabel(char(T(k,8)))
        hold off
        clear D H C
    end
    pause
end
