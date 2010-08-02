% R. McDermott and C. Cruz
% 6-03-2009
% define_drow_variables.m
%
% This script captures the row variables in the 'd line' and stores them in
% a name linked to the header.  The purpose of adding this script is to
% keep 'read_dline.m' as clean as possible.  Note that the cell arrays
% 'parameters' and 'headers' are defined in 'read_dline.m'.
%
% vdir = '../../../Verification/'; for example, is also set in read_dline.

Quantity          = parameters(find(strcmp(headers,'Quantity')));
Group_Key_Label   = parameters(find(strcmp(headers,'Group_Key_Label')));
Group_Style       = parameters(find(strcmp(headers,'Group_Style')));
Fill_Color        = parameters(find(strcmp(headers,'Fill_Color')));
d1_Filename       = [vdir,char(parameters(find(strcmp(headers,'d1_Filename'))))];
d1_Col_Name_Row   = str2num(char(parameters(find(strcmp(headers,'d1_Col_Name_Row')))));
d1_Ind_Col_Name   = char(parameters(find(strcmp(headers,'d1_Ind_Col_Name'))));
d1_Dep_Col_Name   = char(parameters(find(strcmp(headers,'d1_Dep_Col_Name'))));
d1_Data_Row       = str2num(char(parameters(find(strcmp(headers,'d1_Data_Row')))));
d1_Key            = char(parameters(find(strcmp(headers,'d1_Key'))));
d1_Style          = char(parameters(find(strcmp(headers,'d1_Style'))));
d1_Start          = str2num(char(parameters(find(strcmp(headers,'d1_Start')))));
d1_End            = str2num(char(parameters(find(strcmp(headers,'d1_End')))));
d1_Comp_Start     = str2num(char(parameters(find(strcmp(headers,'d1_Comp_Start')))));
d1_Comp_End       = str2num(char(parameters(find(strcmp(headers,'d1_Comp_End')))));
d1_Initial_Value  = str2num(char(parameters(find(strcmp(headers,'d1_Initial_Value')))));
d2_Filename       = [vdir,char(parameters(find(strcmp(headers,'d2_Filename'))))];
d2_Col_Name_Row   = str2num(char(parameters(find(strcmp(headers,'d2_Col_Name_Row')))));
d2_Ind_Col_Name   = char(parameters(find(strcmp(headers,'d2_Ind_Col_Name'))));
d2_Dep_Col_Name   = char(parameters(find(strcmp(headers,'d2_Dep_Col_Name'))));
d2_Data_Row       = str2num(char(parameters(find(strcmp(headers,'d2_Data_Row')))));
d2_Key            = char(parameters(find(strcmp(headers,'d2_Key'))));
d2_Style          = char(parameters(find(strcmp(headers,'d2_Style'))));
d2_Start          = str2num(char(parameters(find(strcmp(headers,'d2_Start')))));
d2_End            = str2num(char(parameters(find(strcmp(headers,'d2_End')))));
d2_Comp_Start     = str2num(char(parameters(find(strcmp(headers,'d2_Comp_Start')))));
d2_Comp_End       = str2num(char(parameters(find(strcmp(headers,'d2_Comp_End')))));
d2_Initial_Value  = str2num(char(parameters(find(strcmp(headers,'d2_Initial_Value')))));
Metric            = char(parameters(find(strcmp(headers,'Metric'))));
Plot_Title        = char(parameters(find(strcmp(headers,'Plot_Title'))));
Ind_Title         = char(parameters(find(strcmp(headers,'Ind_Title'))));
Dep_Title         = char(parameters(find(strcmp(headers,'Dep_Title'))));
Min_Ind           = str2num(char(parameters(find(strcmp(headers,'Min_Ind')))));
Max_Ind           = str2num(char(parameters(find(strcmp(headers,'Max_Ind')))));
Scale_Ind         = str2num(char(parameters(find(strcmp(headers,'Scale_Ind')))));
Min_Dep           = str2num(char(parameters(find(strcmp(headers,'Min_Dep')))));
Max_Dep           = str2num(char(parameters(find(strcmp(headers,'Max_Dep')))));
Scale_Dep         = str2num(char(parameters(find(strcmp(headers,'Scale_Dep')))));
Flip_Axis         = char(parameters(find(strcmp(headers,'Flip_Axis'))));
Title_Position    = str2num(char(parameters(find(strcmp(headers,'Title_Position')))));
Key_Position      = char(parameters(find(strcmp(headers,'Key_Position'))));
Plot_Type         = char(parameters(find(strcmp(headers,'Plot_Type'))));
Plot_Filename     = char(parameters(find(strcmp(headers,'Plot_Filename'))));

