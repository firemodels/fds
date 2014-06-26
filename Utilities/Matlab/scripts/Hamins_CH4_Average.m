%Hamins_CH4.m
%A rewrite of Kevin McGrattan's Hamins_CH4.f90 code into Matlab script.

% Purpose
% =======
% This program computes the time averaged radiative profiles
% Utility program for the FDS validation case Hamins_CH4
%
% This program reads in the devc files for cases Hamins_CH4 1,5,23,21,7,19
% and computes time averaged radial and vertical profiles of radiative flux.
% The averaged proviles are written in files Hamins_CH4_#_devc_aver.csv
%

close all
clear all

% 10 cm burners
mintime = 3.0;
ncol = 37;
nR = 6;
nV = 12;
R = [7.8, 10.0, 12.8, 15.0, 17.5, 20.6];
Z = [0.2, 1.1, 2.2, 3.3, 4.4, 6.6, 8.8, 12.2, 15.2, 20.3, 30.0, 39.3];

% Hamins_CH4_1
infile = 'Hamins_CH4_01_devc.csv';
outfile = 'Hamins_CH4_01_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

% Hamins_CH4_5
infile = 'Hamins_CH4_05_devc.csv';
outfile = 'Hamins_CH4_05_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

% 38 cm burners
mintime = 2.5;
ncol = 41;
nR = 7;
nV = 13;
R = [24.0, 31.2, 40.0, 51.2, 60.0, 70.0, 82.4];
Z = [0.008, 0.044, 0.088, 0.132, 0.176, 0.264, 0.352, 0.488, 0.608, 0.808, 1.200, 1.572, 1.890];
Z = Z*100.0;

% Hamins_CH4_21
infile = 'Hamins_CH4_21_devc.csv';
outfile = 'Hamins_CH4_21_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

% Hamins_CH4_23
infile = 'Hamins_CH4_23_devc.csv';
outfile = 'Hamins_CH4_23_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

% 100 cm burners
mintime = 3.0;
ncol = 41;
nR = 8;
nV = 12;
R = [0.55, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.45];
Z = [0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5];
R = R*100.0;
Z = Z*100.0;

% Hamins_CH4_7
infile = 'Hamins_CH4_07_devc.csv';
outfile = 'Hamins_CH4_07_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

% Hamins_CH4_19
infile = 'Hamins_CH4_19_devc.csv';
outfile = 'Hamins_CH4_19_devc_avg.csv';
Hamins_ReadCH4(infile,outfile,mintime,ncol,nR,nV,R,Z);

cd ..
