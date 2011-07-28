% luetaan spektri .y ASCII tiedostosta
function [a,v]=readyspec(filename)
    temp=importdata(filename); % in L /molcm
    vstart=temp(1)*10^2; % oletetaan että v:n yksiköt 1/cm ja 1/cm->1/m
    vend=temp(2)*10^2;
    npt=temp(3);
    v=linspace(vstart,vend,npt);
    a=temp(4:end)*0.1; % L/mol/cm -> m^3/mol/m
end