% K. McGrattan
% 7-16-2009
% plot_style.m
%
% Preferred style for FDS and Smokeview plots
%
% Usage: When creating a new figure use the following:
%
% plot_style
% figure
% set(gca,'Units',Plot_Units)
% set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
%
% ==>> plot(X,Y,etc)
% ==>> xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
% ==>> ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
%
% set(gca,'FontName',Font_Name)
% set(gca,'FontSize',Label_Font_Size)
%
% ==>> lh=legend(...);
% ==>> set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
%
% set(gcf,'Visible',Figure_Visibility);
% set(gcf,'Units',Paper_Units);
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf',plotname);

style = 'fds'; % set to 'fds', 'paper', etc., as needed (default 'fds')

switch style

    case 'fds'

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
        Scat_Plot_X          = 1.00;
        Scat_Plot_Y          = 0.75;
        Subtitle_Text_Offset = 0.05;

        VerStr_Scale_X   = 0.60;
        VerStr_Scale_Y   = 1.05;

        % Paper properties
        Paper_Units     = 'inches';
        Paper_Width     = 6.5;
        Paper_Height    = 4.5;
        Scat_Paper_Height = 6.0;
        Scat_Paper_Width  = 6.0;

        % Print properties
        Figure_Visibility = 'off';
        Image_File_Type = '-dpdf';

    case 'paper'

        % Font properties
        Font_Name = 'Helvetica'; %get(gca,'fontname')
        Font_Interpreter = 'TeX';
        Key_Font_Size   = 16;
        Title_Font_Size = 20;
        Label_Font_Size = 20;
        Marker_Size = 10;
        Scat_Title_Font_Size = 14;
        Scat_Label_Font_Size = 14;

        % Line properties
        Line_Width      = 1; % get(gca,'linewidth')

        % Plot properties
        Plot_Units      = 'normalized'; %get(gca,'units')
        Pos             = [0.1500    0.1500    0.7750    0.8150]; %get(gca,'position') % [left bottom width height]
        YLabel_Offset   = -0.05; % normalized units
        XLabel_Offset   = -0.05; % normalized units
        Plot_X          = Pos(1);
        Plot_Y          = Pos(2);
        Plot_Width      = Pos(3);
        Plot_Height     = Pos(4); %*.95; % use for exponential notation on y-axis tick labels

        Scat_Plot_Width      = 4.75;
        Scat_Plot_Height     = 4.75;
        Scat_Plot_X          = 0.75;
        Scat_Plot_Y          = 0.75;
        Subtitle_Text_Offset = 0.05;

        VerStr_Scale_X   = 0.60;
        VerStr_Scale_Y   = 1.05;

        % Paper properties
        Paper_Units     = 'inches'; %get(gcf,'paperunits')
        Paper_Pos       = [0.2500    2.5000    8.0000    6.0000]; %get(gcf,'paperposition')
        Paper_Width     = Paper_Pos(3);
        Paper_Height    = Paper_Pos(4);
        Scat_Paper_Height = 6.0;
        Scat_Paper_Width  = 6.0;

        % Print properties
        Figure_Visibility = 'on';
        Image_File_Type = '-dpdf';

    otherwise

end


