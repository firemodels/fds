
% Open directory for plot output
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Define standard plotting parameters
plot_style

% Read the output data.
data = csvread('../../Verification/WUI/vege_mass_conservation_devc.csv', 2)

% The sum of the mass of fuel gas, water vapor, and remaining
% solid material should sum to 1.0 at each time.
mass_delta = 1.0 - sum( data(:,[2:4]), 2 )

plot( data(:,1), mass_delta )
