% K. McGrattan
% 7-16-2009
% plot_style.m
%
% Preferred style for FDS and Smokeview plots 
%
% Some things you might want to set in your script...
%
% set(gca,'Units',Plot_Units)
% set(gca,'FontName',Font_Name)
% set(gca,'FontSize',Title_Font_Size)
% set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
% set(gcf,'DefaultLineLineWidth',Line_Width)

% font properties
Font_Name       = 'Times';
Key_Font_Size   = 12;
Title_Font_Size = 14;
Label_Font_Size = 14;

% line properties
Line_Width      = 1.5;

% plot properties
Plot_Units      = 'inches';
Plot_Width      = 4.5;
Plot_Height     = 3.15;
Plot_X          = 0.75;
Plot_Y          = 0.75;

% paper properties
Paper_Units     = 'inches';
Paper_Width     = 6.0;
Paper_Height    = 4.5;

% print properties
Figure_Visibility = 'on';
