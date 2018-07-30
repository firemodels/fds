% McGrattan
% 7-30-2018
% FM_Burner.m
%
% Read and process FDS output files for FM_Burner cases

close all
clear all

outdir = '../../../out/FM_Burner/FDS_Output_Files/';

fuel_name    = {'C2H4','C3H6','C3H8','CH4'};
res_name     = {'1p0cm'};

for i_fuel=1:4;

   for fds_resolution=1:1;

      DEV = importdata([outdir,'FM_15cm_Burner_',fuel_name{i_fuel},'_devc.csv'],',',2);
      HRR = importdata([outdir,'FM_15cm_Burner_',fuel_name{i_fuel},'_hrr.csv'],',',2);

      Time_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'Time')));
      XO2_FDS  = DEV.data(:,find(strcmp(DEV.colheaders,'"XO2"')));
      Qdot_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'HRR')))-1.;
      Qrad_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'Q_RADI')));
      ntp      = length(Time_FDS);

      fid = fopen([outdir,'FM_15cm_Burner_',fuel_name{i_fuel},'_',res_name{fds_resolution},'.csv'],'wt','n');
      fprintf(fid,'%s\n','XO2,eta,Chi_R');
      for ii=1:ntp
         fprintf(fid,'%5.3f,%6.2f,%6.2f\n',XO2_FDS(ii),Qdot_FDS(ii)/10.,max(0,-Qrad_FDS(ii)/max(0.001,Qdot_FDS(ii))));
      end
      fclose(fid);

   end

end

