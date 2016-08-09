% McDermott
% 8-9-2016
% ht3d_slab.m
%
% Analytical solution to heat transfer in a semi-infinite slab with convection.
% See D. Drysdale, "An Introduction to Fire Dynamcis", 2nd Ed., Wiley, p. 43.
%
% This case uses basically the same parameters as the "convective_cooling"
% solid phase verification case.

close all
clear all

t_end = 1800;
a = .001;
h = 1;
k = 1;
T0 = 1000;
Tinf = 0;

T = @(x,t) T0 + (Tinf-T0) * ( erfc(x./sqrt(a*t)/2) - exp((x*h)/k + a*t/(k/h)^2).*erfc(x./sqrt(a*t)/2 + sqrt(a*t)/(k/h)) );
Ts = @(t) T0 + (Tinf-T0) * ( 1 - exp(a*t/(k/h)^2).*erfc(sqrt(a*t)/(k/h)) );

t_range = linspace(0,t_end,50);

plot(t_range,Ts(t_range),'bo'); hold on

% read FDS results
ddir = '/Volumes/rmcdermo/GitHub/fds-smv_rmcdermo/Verification/Heat_Transfer/';

M = importdata([ddir,'ht3d_slab_devc.csv'],',',2);
jTime = find(strcmp(M.colheaders,'Time'));
jTS = find(strcmp(M.colheaders,'"TS"'));
t_fds = M.data(:,jTime);
Ts_fds = M.data(:,jTS);

plot(t_fds,Ts_fds,'b-')

% % write out analytical ramp at x = 0.5 m

% A = [t_range',T(.5,t_range)'];
% fid = fopen('Tsoln.txt','wt');
% for i=1:length(t_range)
%     fprintf(fid,'%s %f %s %f %s\n','&RAMP ID=''T_ramp'', T=',A(i,1),', F=',A(i,2),' /');
% end
% fclose(fid);

% % write out analytical solution for dataplot

% A = [t_range',Ts(t_range)'];
% fid = fopen([ddir,'ht3d_slab_soln.csv'],'wt');
% fprintf(fid,'%s,%s\n','Time','Ts');
% fprintf(fid,'%s,%s\n','(s)','(C)');
% for i=1:length(t_range)
%     fprintf(fid,'%f,%f\n',A(i,:));
% end
% fclose(fid);
