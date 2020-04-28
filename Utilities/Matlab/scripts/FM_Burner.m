% McGrattan
% 7-30-2018
% FM_Burner.m
%
% Read and process FDS output files for FM_Burner cases

close all
clear all

outdir = '../../../out/FM_Burner/';

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


clear all

outdir = '../../../out/FM_Burner/';
O2_name  = {'20p9','19p0','16p8','15p2'};
res_name = {'2cm','1cm','5mm'};

for i_O2=1:4;

   for fds_resolution=1:3;

      HRR = importdata([outdir,'FM_15cm_Burner_C2H4_',O2_name{i_O2},'_',res_name{fds_resolution},'_hrr.csv'],',',2);

      Time_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'Time')));
      Qdot_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'HRR')));
      Qrad_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'Q_RADI')));
      ntp      = length(Time_FDS);

      chi_r = max(0,-Qrad_FDS./max(0.001,Qdot_FDS));

      fid = fopen([outdir,'FM_15cm_Burner_C2H4_',O2_name{i_O2},'_',res_name{fds_resolution},'_chir.csv'],'wt','n');
      fprintf(fid,'%s\n','Time,Chi_r');
      for ii=2:ntp
         fprintf(fid,'%5.3f,%6.3f\n',Time_FDS(ii),chi_r(ii));
      end
      fclose(fid);

   end
end


expdir = '../../../exp/Submodules/macfp-db/Extinction/FM_Burner/Experimental_Data/Radiation_global/';

CHIR = importdata([expdir,'Radiant_fraction.csv'],',',1);
rfid = fopen([expdir,'Radiant_fraction2.csv'],'wt','n');
fprintf(rfid,'%s\n','Time,20p9 percent,19p0 percent,17p0 percent,15p0 percent');
fprintf(rfid,'%5.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',   0.,CHIR.data(1:4,2));
fprintf(rfid,'%5.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',1000.,CHIR.data(1:4,2));
fclose(rfid);
