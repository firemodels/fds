% McDermott
% 7-20-2009
% turb_model.m
%
% This script generates all plots necessary for the 'Decaying Isotropic
% Turbulence'  section of the verification guide.

energy_decay('csmag_32',32)
energy_decay('csmag_64',64)
energy_decay('dsmag_32',32)
energy_decay('dsmag_64',64)

plotspec('csmag_32',32)
plotspec('csmag_64',64)
plotspec('dsmag_32',32)
plotspec('dsmag_64',64)