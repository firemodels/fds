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
energy_decay('deardorff_32',32)
energy_decay('deardorff_64',64)
energy_decay('vreman_32',32)
energy_decay('vreman_64',64)
energy_decay('rng_32',32)
energy_decay('rng_64',64)
energy_decay('wale_32',32)
energy_decay('wale_64',64)

plotspec('csmag_32',32)
plotspec('csmag_64',64)
plotspec('dsmag_32',32)
plotspec('dsmag_64',64)
plotspec('deardorff_32',32)
plotspec('deardorff_64',64)
plotspec('vreman_32',32)
plotspec('vreman_64',64)
plotspec('rng_32',32)
plotspec('rng_64',64)
plotspec('wale_32',32)
plotspec('wale_64',64)