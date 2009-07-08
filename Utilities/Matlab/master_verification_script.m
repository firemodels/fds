% McDermott
% 5-28-2009
% master_verification_script.m
%
% If you author a section in the verification guide, create a script
% that generates all the graphics (in pdf format) in that section and 
% store the script in the 'scripts' directory.  For example, wall_model.m
% creates all the pdfs for the section on the Werner and Wengle wall
% model.  Also, add your script to the master list below.
%
% If you create and use any customized functions, please store these in
% the 'functions' directory.
%
% To remain backward compatible with the PyroGraph script we have
% included the script read_dline.m.  To utilize this script, add the
% appropriate parameters to a 'd' line in
% verification_data_config_matlab.csv.

cd functions
[saved_data,drange] = dataplot('verification',[2:100]);
cd ..
run scripts/wall_model