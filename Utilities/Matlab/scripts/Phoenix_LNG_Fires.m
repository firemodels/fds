% McGrattan
% 8/19/2021
% Phoenix_LNG_Fires.m
%
% Reformat heat flux output

close all
clear all

outdir = '../../../out/Phoenix_LNG_Fires/';

DEV = importdata([outdir,'Phoenix01_devc.csv'],',',2);
fid = fopen([outdir,'Phoenix01_HF.csv'],'wt','n');
fprintf(fid,'%s\n','Pos-North-Wide,North-Wide,Pos-East-Wide,East-Wide,Pos-South-Wide,South-Wide,Pos-West-Wide,West-Wide');
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',128.6,DEV.data(13,2),108.2,DEV.data(13,5), 92.0,DEV.data(13,8) ,112.4,DEV.data(13,11));
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',178.5,DEV.data(13,3),157.4,DEV.data(13,6),141.2,DEV.data(13,9) ,162.3,DEV.data(13,12));
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',228.3,DEV.data(13,4),207.3,DEV.data(13,7),191.1,DEV.data(13,10),212.1,DEV.data(13,13));
fclose(fid);

clear DEV fid

DEV = importdata([outdir,'Phoenix02_devc.csv'],',',2);
fid = fopen([outdir,'Phoenix02_HF.csv'],'wt','n');
fprintf(fid,'%s\n','Pos-North-Wide,North-Wide,Pos-East-Wide,East-Wide,Pos-South-Wide,South-Wide,Pos-West-Wide,West-Wide');
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',133.0,DEV.data(13,2),103.4,DEV.data(13,5), 87.6,DEV.data(13,8) ,117.7,DEV.data(13,11));
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',182.9,DEV.data(13,3),152.0,DEV.data(13,6),136.8,DEV.data(13,9) ,167.5,DEV.data(13,12));
fprintf(fid,'%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',232.7,DEV.data(13,4),202.1,DEV.data(13,7),186.7,DEV.data(13,10),217.5,DEV.data(13,13));
fclose(fid);


