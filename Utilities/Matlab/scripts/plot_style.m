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
Font_Interpreter = 'LaTeX';
Key_Font_Size   = 12;
Title_Font_Size = 16;
Label_Font_Size = 16;
Scat_Title_Font_Size = 14;
Scat_Label_Font_Size = 14;

% line properties
Line_Width      = 1.0;

% plot properties
Plot_Units      = 'inches';
Plot_Width      = 5.0;
Plot_Height     = 3.4;
Plot_X          = 1.2;
Plot_Y          = 0.8;

Scat_Plot_Width      = 4.75;
Scat_Plot_Height     = 4.75;
Scat_Plot_X          = 0.75;
Scat_Plot_Y          = 0.75;

% paper properties
Paper_Units     = 'inches';
Paper_Width     = 6.5;
Paper_Height    = 4.5;
Scat_Paper_Height = 6.0;
Scat_Paper_Width  = 6.0;

% print properties
Figure_Visibility = 'on';

% svn text position
SVN_Scale_X = 0.80;
SVN_Scale_Y = 1.05;
