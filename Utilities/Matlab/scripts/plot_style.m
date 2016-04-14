% K. McGrattan
% 7-16-2009
% plot_style.m
%
% Preferred style for FDS and Smokeview plots
%
% Usage: When creating a new figure use the following:
%
% figure
% plot_style
% set(gca,'FontName',Font_Name)
% set(gca,'FontSize',Label_Font_Size)
% set(gcf,'Visible',Figure_Visibility);
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% Font properties
Font_Name = 'Times';
Font_Interpreter = 'TeX';
Key_Font_Size   = 12;
Title_Font_Size = 16;
Label_Font_Size = 16;
Scat_Title_Font_Size = 14;
Scat_Label_Font_Size = 14;
Marker_Size = 4;

% Line properties
Line_Width      = 1.0;

% Plot properties
Plot_Units      = 'inches';
Plot_Width      = 5.0;
Plot_Height     = 3.4;
Plot_X          = 1.2;
Plot_Y          = 0.8;

Scat_Plot_Width      = 4.75;
Scat_Plot_Height     = 4.75;
Scat_Plot_X          = 0.75;
Scat_Plot_Y          = 0.75;
Subtitle_Text_Offset = 0.05;

% Paper properties
Paper_Units     = 'inches';
Paper_Width     = 6.5;
Paper_Height    = 4.5;
Scat_Paper_Height = 6.0;
Scat_Paper_Width  = 6.0;

% Print properties
Figure_Visibility = 'on';
Image_File_Type = '-dpdf';
