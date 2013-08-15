% R. McDermott and C. Cruz and S. Hostikka
% 6-06-2012
% dataplot.m
%
% [saved_data, drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir, [drange])
%
% Output:
%    saved_data - cell array containing data needed in scatplot.m
%
%    drange - data range needed for scatplot.m, must be commensurate with
%    saved_data
%
% Input:
%
%    Dataplot_Inputs_File - base configuration file
%
%    Working_Dir - base input file directory
%
%    Manuals_Dir - base plot directory
%
%    [optional] drange - a vector for the 'd' lines you want to read from the
%    config file.  For example, [2:5,7:8,10,12].
%    drange may also be a text string matching the 'Dataname' column in the
%    configuration file.  For example, [currently] the 'WTC' text string
%    identifies drange = [13:18], but is much easier to remember if you
%    happen to work with this set of validation cases frequently.
%    drange input can be overriden by setting the switch_id in the
%    configuration file to "o" in those that you want to be processed only.
%
% Dependencies:
%    ../scripts/define_drow_variables.m
%    dvcread.m
%    parse.m
%
% Example: From the command line within the Matlab/functions/ directory,
%    type
%
%    >> [saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir, [2:4,6:8]);
%
%    >> [saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir, 'WTC');

function [saved_data,drange] = dataplot(varargin)

if nargin<3||nargin>4; 
    display('Error in argument list')
end
if nargin>=3
    Dataplot_Inputs_File = varargin{1};
    Working_Dir = varargin{2};
    Manuals_Dir = varargin{3};
end

% Read in global plot options
plot_style

set(gcf,'DefaultLineLineWidth',Line_Width)
WPos = get(gcf,'Position');
set(gcf,'Position',[WPos(1) WPos(2) 640,420]);
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

% Read configuration file
A = importdata(Dataplot_Inputs_File);
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

n_plots = length(A);

if nargin==4
    drange = varargin{4};
else
    drange = 2:n_plots;
end

if ~isnumeric(drange)
    dataname_col = strcmp(headers,'Dataname');
    dstring = drange;
    drange_index = 0;
    clear drange
else
    dstring = 'null';
end

% Allocate arrays to hold the data for scatterplots
Save_Measured_Metric = zeros(n_plots,10,10);
Save_Predicted_Metric = zeros(n_plots,10,10);

% Search for "o" lines, to process (o)nly those lines.
otest_true = false;
for i=2:n_plots

    if i>n_plots; break; end
    P = textscan(A{i},'%q','delimiter',',');
    parameters = P{:}';
   
    otest = strcmp(parameters(strcmp(headers,'switch_id')),'o');

    if otest
       if ~otest_true
          otest_true = true;
          drange = [];
          dstring = 'null';
       end
       drange = [drange i];
    end
end
   
% Process the "d" or "o" lines one by one
for i=2:n_plots
    
    if i>length(A); break; end

    P = textscan(A{i},'%q','delimiter',',');
    parameters = P{:}';
    
    % Check for shortname specification instead of numeric drange
    if strcmp(dstring,'null')
        itest = ismember(i,drange);
    else
        itest = strcmp(parameters(dataname_col),dstring);
        if itest
            drange_index = drange_index + 1;
            drange(drange_index) = i;
        end
    end
    
    % Check to see if d line has been activated in configuration file
    dtest = strcmp(parameters(strcmp(headers,'switch_id')),'d');

    % Check to see if o line has been activated in configuration file
    otest = strcmp(parameters(strcmp(headers,'switch_id')),'o');
    
    if itest && (dtest || otest)

        figure
        
        define_drow_variables
        
        % Save for scatter plots
        Q1                      = parse(Quantity);
        Save_Quantity(i,1:length(Q1))        = Q1;
        Save_Group_Style(i)     = Group_Style;
        Save_Fill_Color(i)      = Fill_Color;
        Save_Group_Key_Label(i) = Group_Key_Label;
        Save_Dataname(i)        = Stat_Dataname;
        Save_Plot_Filename(i)   = Stat_Plot_Filename;
        Save_Dep_Title(i)       = Stat_Dep_Title;
        Save_Error_Tolerance(i) = Error_Tolerance;
        Save_Metric_Type(i)     = {Metric};
                
        % Plot the experimental data or analytical solution (d1)
        if ~exist(d1_Filename,'file')
           display(['Error: File ', d1_Filename ', does not exist. Skipping case.'])
           continue
        end
        [H M] = dvcread(d1_Filename,d1_Col_Name_Row);
        R1 = parse(d1_Ind_Col_Name);
        S1 = parse(d1_Dep_Col_Name);
        style = parse(d1_Style);
        % Wrap entire d1 dataplot routine in try loop
        % Skips case upon any Matlab error
        try
            for j=1:length(S1)
                d1_Ind_Col = find(strcmp(H,R1(min(j,length(R1)))));
                d1_Dep_Col = find(strcmp(H,S1(j)));
                clear indices
                % Clear flag for stat_x_y metric
                using_stat_x_y = 0;
                using_stat_x_y_check_zero = 0;
                indices = find(d1_Comp_Start<=M(:,d1_Ind_Col) & M(:,d1_Ind_Col)<=d1_Comp_End);
                if strcmp(Metric,'max')
                    Save_Measured_Metric(i,j,1) = max(M(indices,d1_Dep_Col))-d1_Initial_Value;
                elseif strcmp(Metric,'min')
                    Save_Measured_Metric(i,j,1) = d1_Initial_Value-min(M(indices,d1_Dep_Col));
                elseif strcmp(Metric,'maxabs')
                    Save_Measured_Metric(i,j,1) = max(abs(M(indices,d1_Dep_Col)-d1_Initial_Value));
                elseif strfind(Metric,'max_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['max_' '%f' '_' '%f']);
                    if compare_indices(1) == j
                        Save_Measured_Metric(i,1,1) = max(M(indices,d1_Dep_Col))-d1_Initial_Value;
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'mean')
                    Save_Measured_Metric(i,j,1) = abs(mean(M(indices,d1_Dep_Col))-d1_Initial_Value);
                % If mean_x_y is specified for a plot with multiple curves,
                % then get the results from curve x only
                elseif strfind(Metric,'mean_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['mean_' '%f' '_' '%f']);
                    if compare_indices(1) == j
                        Save_Measured_Metric(i,1,1) = abs(mean(M(indices,d1_Dep_Col))-d1_Initial_Value);
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'all')
                    Save_Measured_Metric(i,j,1:length(indices)) = M(indices,d1_Dep_Col)-d1_Initial_Value;
                elseif strcmp(Metric,'threshold')
                    Save_Measured_Metric(i,j,1) = min(M(indices,d1_Dep_Col));
                elseif strcmp(Metric,'area')
                    Save_Measured_Metric(i,j,1) = trapz(M(indices,d1_Ind_Col), M(indices,d1_Dep_Col))-d1_Initial_Value;
                elseif strcmp(Metric,'end')
                    Save_Measured_Metric(i,j,1) = M(indices(end),d1_Dep_Col)-d1_Initial_Value;
                % If end_x_y is specified for a plot with multiple curves,
                % then get the results from curve x only
                elseif strfind(Metric,'end_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['end_' '%f' '_' '%f']);
                    if compare_indices(1) == j
                        Save_Measured_Metric(i,1,1) = M(indices(end),d1_Dep_Col)-d1_Initial_Value;
                        using_stat_x_y_check_zero = 1;
                    end
                else
                    Save_Measured_Metric(i,j,1) = 0;
                end
                % Prevent a value of zero of being returned, which would be erased in statplot using nonzeros()
                if (Save_Measured_Metric(i,j,1) == 0) & (~using_stat_x_y)
                    Save_Measured_Metric(i,j,1) = 1E-12;
                end
                % Special case to pass a zero if using stat_x_y
                if (Save_Measured_Metric(i,1,1) == 0) & (using_stat_x_y_check_zero)
                    Save_Measured_Metric(i,1,1) = 1E-12;
                end
                clear indices
                indices = find(d1_Start<=M(:,d1_Ind_Col) & M(:,d1_Ind_Col)<=d1_End);
                if strcmp(Flip_Axis,'no')
                    X = M(indices,d1_Ind_Col)/Scale_Ind;
                    Y = M(indices,d1_Dep_Col)/Scale_Dep;
                else
                    X = M(indices,d1_Dep_Col)/Scale_Dep;
                    Y = M(indices,d1_Ind_Col)/Scale_Ind;
                end
                if strcmp(Plot_Type,'linear')
                    K(j) = plot(X,Y,char(style(j))); hold on
                elseif strcmp(Plot_Type,'loglog')
                    K(j) = loglog(X,Y,char(style(j))); hold on
                elseif strcmp(Plot_Type,'semilogx')
                    K(j) = semilogx(X,Y,char(style(j))); hold on
                elseif strcmp(Plot_Type,'semilogy')
                    K(j) = semilogy(X,Y,char(style(j))); hold on
                end
            end
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname, '); check syntax of analytical/expected/experimental (d1) columns. Skipping case.'])
            continue
        end
        
        % Plot the FDS or model data (d2)
        if ~exist(d2_Filename,'file')
           display(['Error: File ', d2_Filename, ' does not exist. Skipping case.'])
           continue
        end
        [H M] = dvcread(d2_Filename,d2_Col_Name_Row);
        R2 = parse(d2_Ind_Col_Name);
        S2 = parse(d2_Dep_Col_Name);
        style = parse(d2_Style);
        % Wrap entire d2 dataplot routine in try loop
        % Skips case upon any Matlab error
        try
            for j=1:length(S2)
                d2_Ind_Col = find(strcmp(H,R2(min(j,length(R2)))));
                d2_Dep_Col = find(strcmp(H,S2(j)));
                clear indices
                % Clear flag for stat_x_y metric
                using_stat_x_y = 0;
                using_stat_x_y_check_zero = 0;
                indices = find(d2_Comp_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_Comp_End);
                if strcmp(Metric,'max')
                    Save_Predicted_Metric(i,j,1) = max(M(indices,d2_Dep_Col))-d2_Initial_Value;
                elseif strcmp(Metric,'min')
                    Save_Predicted_Metric(i,j,1) = d2_Initial_Value-min(M(indices,d2_Dep_Col));
                elseif strcmp(Metric,'maxabs')
                    Save_Predicted_Metric(i,j,1) = max(abs(M(indices,d2_Dep_Col)-d2_Initial_Value));
                elseif strfind(Metric,'max_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['max_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = max(M(indices,d2_Dep_Col))-d2_Initial_Value;
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'mean')
                    Save_Predicted_Metric(i,j,1) = abs(mean(M(indices,d2_Dep_Col))-d2_Initial_Value);
                % If mean_x_y is specified for a plot with multiple curves,
                % then get the results from curve y only
                elseif strfind(Metric,'mean_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['mean_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = abs(mean(M(indices,d2_Dep_Col))-d2_Initial_Value);
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'all')
                    Save_Predicted_Metric(i,j,1:length(indices)) = M(indices,d2_Dep_Col)-d2_Initial_Value;
                elseif strcmp(Metric,'threshold')
                    Save_Predicted_Metric(i,j,1) = min(M(indices,d2_Dep_Col))-d2_Initial_Value;
                elseif strcmp(Metric,'area')
                    Save_Predicted_Metric(i,j,1) = trapz(M(indices,d2_Ind_Col), M(indices,d2_Dep_Col))-d2_Initial_Value;
                elseif strcmp(Metric,'end')
                    Save_Predicted_Metric(i,j,1) = M(indices(end),d2_Dep_Col)-d2_Initial_Value;
                % If end_x_y is specified for a plot with multiple curves,
                % then get the results from curve y only
                elseif strfind(Metric,'end_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['end_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = M(indices(end),d2_Dep_Col)-d2_Initial_Value;
                        using_stat_x_y_check_zero = 1;
                    end
                else
                    Save_Predicted_Metric(i,j,1) = 0;
                end
                % Prevent a value of zero of being returned, which would be erased in statplot using nonzeros()
                if (Save_Predicted_Metric(i,j,1) == 0) & (~using_stat_x_y)
                    Save_Predicted_Metric(i,j,1) = 1E-12;
                end
                % Special case to pass a zero if using stat_x_y
                if (Save_Predicted_Metric(i,1,1) == 0) & (using_stat_x_y_check_zero)
                    Save_Predicted_Metric(i,1,1) = 1E-12;
                end
                clear indices
                indices = find(d2_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_End);
                if strcmp(Flip_Axis,'no')
                    X = M(indices,d2_Ind_Col)/Scale_Ind;
                    Y = M(indices,d2_Dep_Col)/Scale_Dep;
                else
                    X = M(indices,d2_Dep_Col)/Scale_Dep;
                    Y = M(indices,d2_Ind_Col)/Scale_Ind;
                end
                if strcmp(Plot_Type,'linear')
                    K(length(S1)+j) = plot(X,Y,char(style(j)));
                elseif strcmp(Plot_Type,'loglog')
                    K(length(S1)+j) = loglog(X,Y,char(style(j)));
                elseif strcmp(Plot_Type,'semilogx')
                    K(length(S1)+j) = semilogx(X,Y,char(style(j)));
                elseif strcmp(Plot_Type,'semilogy')
                    K(length(S1)+j) = semilogy(X,Y,char(style(j)));
                end
            end
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname, '); check syntax of FDS/model results (d2) columns. Skipping case.'])
            continue
        end
        hold off
        
        % Wrap entire plot/save routine in try loop
        % Skips case upon any Matlab error
        try
            if strcmp(Plot_Type,'linear') & strcmp(Flip_Axis,'no')
                X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
                Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
            elseif strcmp(Plot_Type,'linear') & strcmp(Flip_Axis,'yes')
                X_Title_Position = Min_Dep+Title_Position(1)*(Max_Dep-Min_Dep);
                Y_Title_Position = Min_Ind+Title_Position(2)*(Max_Ind-Min_Ind);
            elseif strcmp(Plot_Type,'loglog') & strcmp(Flip_Axis,'no')
                X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
                Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));
            elseif strcmp(Plot_Type,'loglog') & strcmp(Flip_Axis,'yes')
                X_Title_Position = 10^(log10(Min_Dep)+Title_Position(1)*(log10(Max_Dep)-log10(Min_Dep)));
                Y_Title_Position = 10^(log10(Min_Ind)+Title_Position(2)*(log10(Max_Ind)-log10(Min_Ind)));
            elseif strcmp(Plot_Type,'semilogx') & strcmp(Flip_Axis,'no')
                X_Title_Position = 10^(log10(Min_Ind)+Title_Position(1)*(log10(Max_Ind)-log10(Min_Ind)));
                Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
            elseif strcmp(Plot_Type,'semilogx') & strcmp(Flip_Axis,'yes')
                X_Title_Position = 10^(log10(Min_Dep)+Title_Position(1)*(log10(Max_Dep)-log10(Min_Dep)));
                Y_Title_Position = Min_Dep+Title_Position(2)*(Max_Dep-Min_Dep);
            elseif strcmp(Plot_Type,'semilogy') & strcmp(Flip_Axis,'no')
                X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
                Y_Title_Position = 10^(log10(Min_Dep)+Title_Position(2)*(log10(Max_Dep)-log10(Min_Dep)));
            elseif strcmp(Plot_Type,'semilogy') & strcmp(Flip_Axis,'yes')
                X_Title_Position = Min_Ind+Title_Position(1)*(Max_Ind-Min_Ind);
                Y_Title_Position = 10^(log10(Min_Ind)+Title_Position(2)*(log10(Max_Ind)-log10(Min_Ind)));
            end

            set(gca,'FontName',Font_Name)
            set(gca,'FontSize',Label_Font_Size)

            if strcmp(Flip_Axis,'no')
                xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
                ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
                axis([Min_Ind Max_Ind Min_Dep Max_Dep])
                text(X_Title_Position,Y_Title_Position,...
                    Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
            else
                xlabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
                ylabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
                axis([Min_Dep Max_Dep Min_Ind Max_Ind])
                text(X_Title_Position,Y_Title_Position,...
                    Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
            end
            if size(Key_Position)>0
                legend_handle = legend(K,[parse(d1_Key),parse(d2_Key)],'Location',Key_Position);
                if strcmp(Key_Position,'EastOutside')
                   pos = get(legend_handle,'position');
                   set(legend_handle,'position',[Paper_Width (Plot_Y+(Plot_Height-pos(4))/2) pos(3:4)])
                end
                if strcmp(Key_Position,'SouthEastOutside')
                   pos = get(legend_handle,'position');
                   set(legend_handle,'position',[Paper_Width Plot_Y pos(3:4)])
                end
                set(legend_handle,'Interpreter',Font_Interpreter);
                set(legend_handle,'Fontsize',Key_Font_Size);
                set(legend_handle,'Box','on');
                if size(d1_Tick)>0
                   set(gca,'XTick',d1_Tick)
                end
                if size(d2_Tick)>0
                   set(gca,'YTick',d2_Tick)
                end
                if size(Legend_XYWidthHeight)>0
                   legend_position=get(legend_handle,'Position');
                   if Legend_XYWidthHeight(1)>0; legend_position(1)=Legend_XYWidthHeight(1); end % X
                   if Legend_XYWidthHeight(2)>0; legend_position(2)=Legend_XYWidthHeight(2); end % Y
                   if Legend_XYWidthHeight(3)>0; legend_position(3)=Legend_XYWidthHeight(3); end % Width
                   if Legend_XYWidthHeight(4)>0; legend_position(4)=Legend_XYWidthHeight(4); end % Height
                   set(legend_handle,'Position',legend_position)
                end
            end

            % Add SVN if file is available
            if exist(SVN_Filename,'file')
                SVN = importdata(SVN_Filename);
                x_lim = get(gca,'XLim');
                y_lim = get(gca,'YLim');
                if strcmp(Plot_Type,'loglog')
                    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
                    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
                elseif strcmp(Plot_Type,'semilogx')
                    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
                    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
                elseif strcmp(Plot_Type,'semilogy')
                    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
                    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
                else
                    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
                    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
                end
                text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
                    'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
            end

            % Save plot file
            PDF_Paper_Width = Paper_Width_Factor*Paper_Width;

            set(gcf,'Visible',Figure_Visibility);
            set(gcf,'PaperUnits',Paper_Units);
            set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
            set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Paper_Height]); 
            display(['Printing plot ',num2str(i),'...'])
            print(gcf,Image_File_Type,[Manuals_Dir,Plot_Filename])
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname, '); check syntax of plot/save settings. Skipping case.'])
            continue
        end    
        
    end
    clear S1 S2 K style H M X Y P parameters
    close all
end

clear A

% Pack data for use in scatplot
saved_data = [{Save_Quantity'},...
              {Save_Group_Style'},...
              {Save_Fill_Color'},...
              {Save_Group_Key_Label'},...
              {Save_Measured_Metric},...
              {Save_Predicted_Metric},...
              {Save_Dataname'},...
              {Save_Plot_Filename'},...
              {Save_Dep_Title'},...
              {Save_Error_Tolerance'},...
              {Save_Metric_Type'}];

display('dataplot completed successfully!')


