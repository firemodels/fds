% McGrattan
% 2-26-2018
% umd_line_burner_process.m
%
% Read and process FDS output files for UMD Line Burner

close all
clear all

Lf_pts = 50;
Lf_min = 0;
Lf_max = 1;
Lf_dz  = (Lf_max-Lf_min)/(Lf_pts-1);
Lf_z   = Lf_min:Lf_dz:Lf_max;
Lf_fac = 0.98;
Lf_dt  = 10;

outdir = '../../../out/UMD_Line_Burner/';

fuel_name    = {'methane','propane'};
res_name     = {'1p25cm','p625cm','p3125cm'};

for i_fuel=1:2;

   for fds_resolution=1:3;

      DEV = importdata([outdir,fuel_name{i_fuel},'_dx_',res_name{fds_resolution},'_devc.csv'],',',2);
      HRR = importdata([outdir,fuel_name{i_fuel},'_dx_',res_name{fds_resolution},'_hrr.csv'],',',2);

      Time_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'Time')));
      XO2_FDS  = DEV.data(:,find(strcmp(DEV.colheaders,'"XO2"')));
      Qdot_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'HRR')));
      Qrad_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'Q_RADI')));
      q_R_FDS  = 0.5*( DEV.data(:,find(strcmp(DEV.colheaders,'"qrad1"'))) + DEV.data(:,find(strcmp(DEV.colheaders,'"qrad2"'))) );
      colLf01  = find(strcmp(DEV.colheaders,'"Lf-01"'));
      colLf50  = find(strcmp(DEV.colheaders,'"Lf-50"'));
      ntp      = length(Time_FDS);
      Lf_FDS   = zeros(1,ntp);

      for n=1:ntp
          hrrpul = DEV.data(n,colLf01:colLf50);
          q_total = sum(hrrpul);
          q_sum = 0.;
          for i=1:Lf_pts
              q_sum = q_sum+hrrpul(i);
              if q_sum>Lf_fac*q_total
                  Lf_FDS(n) = Lf_z(i);
                  break;
              end
          end
      end
      Lf_tmp = Lf_FDS;
      for n=1:ntp
          indx_range = [find(Time_FDS>(Time_FDS(n)-Lf_dt),1):n];
          Lf_FDS(n) = mean(Lf_tmp(indx_range));
      end

      fid = fopen([outdir,fuel_name{i_fuel},'_dx_',res_name{fds_resolution},'.csv'],'wt','n');
      fprintf(fid,'%s\n','XO2,eta,Chi_R,Lf,q_R');
      for ii=1:ntp
         fprintf(fid,'%5.3f,%6.2f,%6.2f,%6.2f,%6.2f\n',XO2_FDS(ii),Qdot_FDS(ii)/50.,max(0,-Qrad_FDS(ii)/max(0.001,Qdot_FDS(ii))),Lf_FDS(ii),q_R_FDS(ii));
      end
      fclose(fid);

   end

end

