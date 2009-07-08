% R. McDermott and C. Cruz
% 6-11-2009
% define_qrow_variables.m
%
% This script captures the row variables in the 'q' line and stores them in
% a name linked to the header.  The purpose of adding this script is to
% keep 'scatterplot.m' as clean as possible.  Note that the cell arrays
% 'parameters' and 'headers' are defined in 'scatterplot.m'.

P = textscan(Q{j},'%q','delimiter',',');
parameters = P{:}';

Quantity_Label      = char(parameters(find(strcmp(headers,'Quantity_Label'))));
Scatter_Plot_Title  = char(parameters(find(strcmp(headers,'Scatter_Plot_Title'))));
Ind_Title           = char(parameters(find(strcmp(headers,'Ind_Title'))));
Dep_Title           = char(parameters(find(strcmp(headers,'Dep_Title'))));
Plot_Min            = str2num(char(parameters(find(strcmp(headers,'Plot_Min')))));
Plot_Max            = str2num(char(parameters(find(strcmp(headers,'Plot_Max')))));
Title_Position      = str2num(char(parameters(find(strcmp(headers,'Title_Position')))));
Key_Position        = char(parameters(find(strcmp(headers,'Key_Position'))));
Plot_Width          = str2num(char(parameters(find(strcmp(headers,'Plot_Width')))));
Sigma_2_E           = str2num(char(parameters(find(strcmp(headers,'Sigma_2_E')))));
Model_Error         = char(parameters(find(strcmp(headers,'Model_Error'))));
Plot_Filename       = char(parameters(find(strcmp(headers,'Plot_Filename'))));