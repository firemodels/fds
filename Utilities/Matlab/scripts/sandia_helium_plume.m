% McDermott
% 2-15-10
% sandia_helium_plume.m
%
% Plot Sandia 1m helium plume data and FDS results

close all
clear all

addpath '..\..\Validation\Sandia_He_1m\Experimental_Data'
addpath '..\..\Validation\Sandia_He_1m\FDS_Input_Files'

radial_profile('W2.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_W2',[2:22],0.2,'bo','b-','vertical velocity (m/s)',-.5,.5,.1,0,5,1)
radial_profile('W4.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_W4',[23:43],0.2,'bo','b-','vertical velocity (m/s)',-.5,.5,.1,0,5,1)
radial_profile('W6.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_W6',[44:64],0.2,'bo','b-','vertical velocity (m/s)',-.5,.5,.1,0,5,1)

radial_profile('U2.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_U2',[65:85],0.4,'ro','r-','horizontal velocity (m/s)',-.5,.5,.1,-1,1,.25)
radial_profile('U4.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_U4',[86:106],0.4,'ro','r-','horizontal velocity (m/s)',-.5,.5,.1,-1,1,.25)
radial_profile('U6.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_U6',[107:127],0.4,'ro','r-','horizontal velocity (m/s)',-.5,.5,.1,-1,1,.25)

radial_profile('YHe2.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_YHe2',[128:148],0.4,'go','g-','helium mass fraction',-.5,.5,.1,0,1,.2)
radial_profile('YHe4.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_YHe4',[149:169],0.4,'go','g-','helium mass fraction',-.5,.5,.1,0,1,.2)
radial_profile('YHe6.csv','Sandia_He_1m_devc.csv','Sandia_He_1m\Sandia_He_1m_YHe6',[170:190],0.4,'go','g-','helium mass fraction',-.5,.5,.1,0,1,.2)