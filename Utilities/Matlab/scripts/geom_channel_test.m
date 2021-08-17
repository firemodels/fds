% McDermott
% 8-11-2021
% geom_channel_test.m
%
% Temporary home for checking AREA integration for geom_channel.fds

close all
clear all

M = importdata('../../Verification/Complex_Geometry/geom_channel_devc.csv',',',2);

U0 = M.data(end,find(strcmp(M.colheaders,'"U0"')));
A0 = M.data(end,find(strcmp(M.colheaders,'"A0"')));
U1 = M.data(end,find(strcmp(M.colheaders,'"U1"')));
A1 = M.data(end,find(strcmp(M.colheaders,'"A1"')));
U2 = M.data(end,find(strcmp(M.colheaders,'"U2"')));
A2 = M.data(end,find(strcmp(M.colheaders,'"A2"')));
U3 = M.data(end,find(strcmp(M.colheaders,'"U3"')));
A3 = M.data(end,find(strcmp(M.colheaders,'"A3"')));

% check errors
error_tolerance = 1e-6;

A0_error = abs(A0-1); if A0_error>error_tolerance; display(['Matlab Warning: geom_channel.fds A0_error = ',num2str(A0_error)]); end
A1_error = abs(A1-1); if A1_error>error_tolerance; display(['Matlab Warning: geom_channel.fds A1_error = ',num2str(A1_error)]); end
A2_error = abs(A2-1); if A2_error>error_tolerance; display(['Matlab Warning: geom_channel.fds A2_error = ',num2str(A2_error)]); end
A3_error = abs(A3-1); if A3_error>error_tolerance; display(['Matlab Warning: geom_channel.fds A3_error = ',num2str(A3_error)]); end