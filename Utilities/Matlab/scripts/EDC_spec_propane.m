function dy=EDC_spec_propane(t,y)

% (1) - nitrogen
% (2) - propane
% (3) - oxygen
% (4) - co2
% (5) - h2o

%-----------------------
% Initialize variables
%-----------------------

m=1;                % chemical reaction coefficient fuel
s=5;                % chemical reaction coefficient oxygen
c=3;                % chemical reaction coefficient carbon dioxide
h=4;                % chemical reaction coefficient water vapor
tau_mix=0.125;      % set mixing time [s]

%-----------------------
% Logic for setting molar production rate
% No suppression in this case
%-----------------------
r1=min([y(2) y(3)/s]);
r_co2(1)=m*r1/tau_mix;

%-----------------------
% Integration
%-----------------------
dy(1)= 0.0;
dy(2)=-(1/m)*r_co2(1);
dy(3)=-(s/m)*r_co2(1);
dy(4)= c*r_co2(1);
dy(5)= h*r_co2(1);

dy=[dy(1);dy(2);dy(3);dy(4);dy(5)];