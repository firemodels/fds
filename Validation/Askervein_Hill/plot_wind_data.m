% McDermott
% 12-9-09
% plot_wind_data.m

close all
clear all

% measurements
xdevc = [414,651,763,858,920,984,1055,1124,1262];
ydevc = [1080,1336,1456,1562,1625,1695,1770,1842,1975];
r = sign(xdevc-xdevc(6)).*sqrt((xdevc-xdevc(6)).^2 + (ydevc-ydevc(6)).^2); % point 6 is hill top

S1 = xlsread('Table1.3_A_data.xls','d9:d17');
S2 = xlsread('Table1.3_A_data.xls','d23:d31');

H(1)=plot(r,S1,'b-'); hold on
plot(r,S2,'b-')

% FDS results
M = csvread('Askervein_devc.csv',2,0);
range = find(M(:,1)>0); %length(M(:,1));
ASW85 = mean(M(range,2));
ASW50 = mean(M(range,3));
ASW35 = mean(M(range,4));
ASW20 = mean(M(range,5));
ASW10 = mean(M(range,6));
HT    = mean(M(range,7));
ANE10 = mean(M(range,8));
ANE20 = mean(M(range,9));
ANE40 = mean(M(range,10));
FDS = [ASW85,ASW50,ASW35,ASW20,ASW10,HT,ANE10,ANE20,ANE40];
H(2)=plot(r,FDS,'ro-');
xlabel('distance from hill top (m)')
ylabel('wind speed (m/s)')
legend(H,'measurements','FDS results','location','northwest')
axis([-1000 400 0 20])

% print('-dpdf','A_fds_csmag.pdf')
% csvwrite('A_fds_csmag.csv',[r',FDS'])
% 
% figure
% M = csvread('ASW10_csmag_vel_profile.csv',2,0);
% range = find(M(:,3)>0);
% z_mean = M(range,2);
% s_mean = M(range,3);
% H(2)=plot(s_mean,z_mean,'bo-'); hold on
% M = csvread('ASW10_csmag_vel_instant_profile.csv',2,0);
% range = find(M(:,3)>0);
% z_inst = M(range,2);
% s_inst = M(range,3);
% H(1)=plot(s_inst,z_inst,'r*-');
% xlabel('wind speed (m/s)')
% ylabel('height (m)')
% legend(H,'instantaneous','mean','location','northwest')
% 
% print('-dpdf','ASW10_profiles_csmag.pdf')
% csvwrite('ASW10_profiles_csmag.csv',[z_mean,s_mean,s_inst])

