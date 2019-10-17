% R. McDermott and C. Cruz and S. Hostikka
% 6-06-2012
% dataplot.m
%
% Detailed instructions can be found here:
% https://github.com/firemodels/fds/wiki/Using-the-Matlab-script-dataplot.m
%
% [saved_data, drange] = dataplot(Dataplot_Inputs_File, EXP_Dir, OUT_Dir, Manuals_Dir, [drange])
%
% Output:
%
%    saved_data - cell array containing data needed in scatplot.m
%
%    drange - data range needed for scatplot.m, must be commensurate with
%    saved_data
%
% Input:
%
%    Dataplot_Inputs_File - base configuration file
%
%    EXP_Dir - location of experimental data repository
%
%    OUT_Dir - location of output file repository
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
%
%    ../scripts/define_drow_variables.m
%    dvcread.m
%    parsepipe.m
%    parseplus.m
%
% Example: From the command line within the Matlab/scripts/ directory,
%    type
%
%    >> [saved_data,drange] = dataplot(Dataplot_Inputs_File, EXP_Dir, OUT_Dir, Manuals_Dir, [2:4,6:8]);
%
%    >> [saved_data,drange] = dataplot(Dataplot_Inputs_File, EXP_Dir, OUT_Dir, Manuals_Dir, 'WTC');
%
% Special switch_id tags:
%
%    'd' -- Process this data line as usual (exception: see 'o' below)
%
%    's' -- Skip this line
%
%    'o' -- Add 'o' in the switch_id column (first column) of FDS_validation_dataplot_inputs.csv to process "only" these lines.
%
%    'f' -- Follow the previous line and "hold on" the figure window, adding this line to the current plot.
%
%    'g' -- Generate plot, but ignore in scatplot.  Good for cases under development.

function [saved_data,drange] = dataplot(varargin)

if nargin<4||nargin>5;
    display('Error in argument list')
end
if nargin>=4
    Dataplot_Inputs_File = varargin{1};
    EXP_Dir = varargin{2};
    OUT_Dir = varargin{3};
    Manuals_Dir = varargin{4};
end

% Read in global plot options
plot_style

% Read configuration file
A = importdata(Dataplot_Inputs_File);
H = textscan(A{1},'%q','delimiter',',');
headers = H{:}'; clear H

n_plots = length(A);

if nargin==5
    drange = varargin{5};
else
    drange = 2:n_plots;
end

drange_index = 0;
if ~isnumeric(drange)
    dataname_col = strcmp(headers,'Dataname');
    dstring = drange;
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
    end
    if any(itest)
        drange_index = drange_index + 1;
        drange(drange_index) = i;
    end

    % Check to see if d line has been activated in configuration file
    dtest = strcmp(parameters(strcmp(headers,'switch_id')),'d');

    % Check to see if o line has been activated in configuration file
    otest = strcmp(parameters(strcmp(headers,'switch_id')),'o');

    % Check to see if f line has been activated in configuration file
    ftest = strcmp(parameters(strcmp(headers,'switch_id')),'f'); % used for multiple lines on same plot

    % Check to see if g line has been activated in configuration file
    gtest = strcmp(parameters(strcmp(headers,'switch_id')),'g'); % used to ignore scatplot

    if any(itest) && (dtest || otest || ftest || gtest)

        % remove this plot from drange if gtest
        if (gtest)
            drange(drange_index)=[];
            drange_index=drange_index-1;
        end

        if ~ftest
            if exist('K1')
                clear K1
            end
            if exist('K2')
                clear K2
            end
            if exist('d1_Key')
                clear d1_Key
            end
            if exist('d2_Key')
                clear d2_Key
            end
            close all
            f1 = figure;
            set(f1,'Visible',Figure_Visibility)
        else
            hold on
            try
                K1_save = K1;
                K2_save = K2;
                d1_Key_save = d1_Key;
                d2_Key_save = d2_Key;
            catch
                display(['Error: Problem with dataplot row ', num2str(i),'.  Skipping case.'])
                continue
            end
        end
        set(gca,'Units',Plot_Units)
        set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

        define_drow_variables

        % Save for scatter plots
        Q1                      = parsepipe(Quantity);
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
        [H M] = dvcread(d1_Filename,d1_Col_Name_Row,d1_Data_Row);
        R1 = parsepipe(d1_Ind_Col_Name);
        S1 = parsepipe(d1_Dep_Col_Name);
        style = parsepipe(d1_Style);
        % Wrap entire d1 dataplot routine in try loop
        % Skips case upon any Matlab error
        try
            for j=1:length(S1)
                d1_Ind_Col = find(strcmp(strtrim(H),strtrim(R1(min(j,length(R1))))));
                d1_Dep_Col = find(strcmp(strtrim(H),strtrim(S1(j))));
                Save_Measured_Quantity(i,j) = S1(j);
                clear indices
                % Clear flag for stat_x_y metric
                using_stat_x_y = 0;
                using_stat_x_y_check_zero = 0;
                indices = find(d1_Comp_Start    <=M(:,d1_Ind_Col)    & M(:,d1_Ind_Col)   <=d1_Comp_End & ...
                               d1_Dep_Comp_Start<=M(:,d1_Dep_Col(1)) & M(:,d1_Dep_Col(1))<=d1_Dep_Comp_End);
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
                elseif strcmp(Metric,'slope')
                    p = polyfit(M(indices,d1_Ind_Col),M(indices,d1_Dep_Col),1);
                    Save_Measured_Metric(i,j,1) = p(1);
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
                        Save_Measured_Quantity(i,1) = S1(j);
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
                if ~ftest
                    if strcmp(Plot_Type,'linear')
                        K1(j) = plot(X,Y,char(style(j))); hold on
                        if strcmp(Metric,'slope') plot([0 10000],[p(2),p(2)+10000*p(1)],'r-'); end
                    elseif strcmp(Plot_Type,'loglog')
                        K1(j) = loglog(X,Y,char(style(j))); hold on
                    elseif strcmp(Plot_Type,'semilogx')
                        K1(j) = semilogx(X,Y,char(style(j))); hold on
                    elseif strcmp(Plot_Type,'semilogy')
                        K1(j) = semilogy(X,Y,char(style(j))); hold on
                    end
                    set(K1(j),'linewidth',Line_Width)
                else
                    if ~strcmp(char(style(j)),'blank')
                       if strcmp(Plot_Type,'linear')
                           K1(length(K1_save)+j) = plot(X,Y,char(style(j))); hold on
                           if strcmp(Metric,'slope') plot([0 10000],[p(2),p(2)+10000*p(1)],'r-'); end
                       elseif strcmp(Plot_Type,'loglog')
                           K1(length(K1_save)+j) = loglog(X,Y,char(style(j))); hold on
                       elseif strcmp(Plot_Type,'semilogx')
                           K1(length(K1_save)+j) = semilogx(X,Y,char(style(j))); hold on
                       elseif strcmp(Plot_Type,'semilogy')
                           K1(length(K1_save)+j) = semilogy(X,Y,char(style(j))); hold on
                       end
                       set(K1(length(K1_save)+j),'linewidth',Line_Width)
                    end
                end
            end
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname,...
                '); check syntax of analytical/expected/experimental (d1) columns. Skipping case.'])
            continue
        end

        % Plot the FDS or model data (d2)
        if ~exist(d2_Filename,'file')
           display(['Error: File ', d2_Filename, ' does not exist. Skipping case.'])
           continue
        end
        [H M] = dvcread(d2_Filename,d2_Col_Name_Row,d2_Data_Row);
        R2 = parsepipe(d2_Ind_Col_Name);
        S2 = parsepipe(d2_Dep_Col_Name);
        style = parsepipe(d2_Style);
        % Wrap entire d2 dataplot routine in try loop
        % Skips case upon any Matlab error
        try
            for j=1:length(S2)

                if strcmp(char(style(j)),'none')
                    continue
                end

                % check for "+" operator on columns (see hrrpuv_reac for examples)
                SP = parseplus(S2(j));
                Save_Predicted_Quantity(i,j) = S2(j);
                d2_Ind_Col = find(strcmp(strtrim(H),strtrim(R2(min(j,length(R2))))));
                for jj=1:length(SP)
                    d2_Dep_Col(jj) = find(strcmp(strtrim(H),strtrim(SP(jj))));
                end
                clear indices

                % Clear flag for stat_x_y metric
                using_stat_x_y = 0;
                using_stat_x_y_check_zero = 0;
                indices = find(d2_Comp_Start    <=M(:,d2_Ind_Col)    & M(:,d2_Ind_Col)   <=d2_Comp_End & ...
                               d2_Dep_Comp_Start<=M(:,d2_Dep_Col(1)) & M(:,d2_Dep_Col(1))<=d2_Dep_Comp_End);

                M_Ind = M(indices,d2_Ind_Col);
                M_Dep = sum(M(indices,d2_Dep_Col),2);

                if strcmp(Metric,'max')
                    Save_Predicted_Metric(i,j,1) = max(M_Dep)-d2_Initial_Value;
                elseif strcmp(Metric,'min')
                    Save_Predicted_Metric(i,j,1) = d2_Initial_Value-min(M_Dep);
                elseif strcmp(Metric,'maxabs')
                    Save_Predicted_Metric(i,j,1) = max(abs(M_Dep-d2_Initial_Value));
                elseif strfind(Metric,'max_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['max_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = max(M_Dep)-d2_Initial_Value;
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'slope')
                    p = polyfit(M(indices,d2_Ind_Col),M(indices,d2_Dep_Col),1);
                    Save_Predicted_Metric(i,j,1) = p(1);
                elseif strcmp(Metric,'mean')
                    Save_Predicted_Metric(i,j,1) = abs(mean(M_Dep)-d2_Initial_Value);
                % If mean_x_y is specified for a plot with multiple curves,
                % then get the results from curve y only
                elseif strfind(Metric,'mean_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['mean_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = abs(mean(M_Dep)-d2_Initial_Value);
                        using_stat_x_y_check_zero = 1;
                    end
                elseif strcmp(Metric,'all')
                    Save_Predicted_Metric(i,j,1:length(indices)) = M_Dep-d2_Initial_Value;
                elseif strcmp(Metric,'threshold')
                    Save_Predicted_Metric(i,j,1) = min(M_Dep)-d2_Initial_Value;
                elseif strcmp(Metric,'area')
                    Save_Predicted_Metric(i,j,1) = trapz(M_Ind,M_Dep)-d2_Initial_Value;
                elseif strcmp(Metric,'end')
                    Save_Predicted_Metric(i,j,1) = M_Dep(end)-d2_Initial_Value;
                % If end_x_y is specified for a plot with multiple curves,
                % then get the results from curve y only
                elseif strfind(Metric,'end_')
                    using_stat_x_y = 1;
                    compare_indices = sscanf(Metric, ['end_' '%f' '_' '%f']);
                    if compare_indices(2) == j
                        Save_Predicted_Metric(i,1,1) = M_Dep(end)-d2_Initial_Value;
                        Save_Predicted_Quantity(i,1) = S2(j);
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

                % Plots
                indices = find(d2_Start<=M(:,d2_Ind_Col) & M(:,d2_Ind_Col)<=d2_End);
                M_Ind = M(indices,d2_Ind_Col);
                M_Dep = sum(M(indices,d2_Dep_Col),2);
                clear d2_Dep_Col;
                if strcmp(Flip_Axis,'no')
                    X = M_Ind/Scale_Ind;
                    Y = M_Dep/Scale_Dep;
                else
                    X = M_Dep/Scale_Dep;
                    Y = M_Ind/Scale_Ind;
                end
                if ~ftest
                    if strcmp(Plot_Type,'linear')
                        K2(j) = plot(X,Y,char(style(j)));
                        if strcmp(Metric,'slope') plot([0 10000],[p(2),p(2)+10000*p(1)],'r--'); end
                    elseif strcmp(Plot_Type,'loglog')
                        K2(j) = loglog(X,Y,char(style(j)));
                    elseif strcmp(Plot_Type,'semilogx')
                        K2(j) = semilogx(X,Y,char(style(j)));
                    elseif strcmp(Plot_Type,'semilogy')
                        K2(j) = semilogy(X,Y,char(style(j)));
                    end
                    set(K2(j),'linewidth',Line_Width)
                else
                    if ~strcmp(char(style(j)),'blank')
                       if strcmp(Plot_Type,'linear')
                           K2(length(K2_save)+j) = plot(X,Y,char(style(j)));
                           if strcmp(Metric,'slope') plot([0 10000],[p(2),p(2)+10000*p(1)],'r--'); end
                       elseif strcmp(Plot_Type,'loglog')
                           K2(length(K2_save)+j) = loglog(X,Y,char(style(j)));
                       elseif strcmp(Plot_Type,'semilogx')
                           K2(length(K2_save)+j) = semilogx(X,Y,char(style(j)));
                       elseif strcmp(Plot_Type,'semilogy')
                           K2(length(K2_save)+j) = semilogy(X,Y,char(style(j)));
                       end
                       set(K2(length(K2_save)+j),'linewidth',Line_Width)
                    end
                end
            end
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname,...
                '); check syntax of FDS/model results (d2) columns. Skipping case.'])
            continue
        end

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

            % Inserts title, skips if 'f' switch (avoids overplotting)
          % if ~ftest
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
          % end

            if size(Key_Position)>0
                if ~ftest
                    legend_handle = legend([K1 K2],[parsepipe(d1_Key),parsepipe(d2_Key)],'Location',Key_Position);
                else
                    % this allows us to handle multiple lines on the same plot
                    if ~strcmp(d1_Key,'blank')
                       d1_Key = [d1_Key_save,'|',d1_Key];
                    else
                       d1_Key = d1_Key_save;
                    end
                    if ~strcmp(d2_Key,'blank')
                       d2_Key = [d2_Key_save,'|',d2_Key];
                    else
                       d2_Key = d2_Key_save;
                    end
                    legend_handle = legend([K1 K2],[parsepipe(d1_Key),parsepipe(d2_Key)],'Location',Key_Position);
                end
                if strcmp(Key_Position,'EastOutside')
                   set(legend_handle,'Units',Paper_Units)
                   pos = get(legend_handle,'position');
                   set(legend_handle,'position',[Paper_Width pos(2:4)])
                end
                if strcmp(Key_Position,'SouthEastOutside')
                   set(legend_handle,'Units',Paper_Units)
                   pos = get(legend_handle,'position');
                   set(legend_handle,'position',[Paper_Width pos(2:4)])
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

            % Add version string if file is available, skips if 'f' switch (avoids overplotting)
            if ~ftest
                addverstr(gca,VerStr_Filename,Plot_Type,VerStr_Scale_X,VerStr_Scale_Y,Font_Name,Font_Interpreter)
            end

            % Save plot file
            PDF_Paper_Width = Paper_Width_Factor*Paper_Width;

            set(gcf,'Visible',Figure_Visibility);
            set(gcf,'PaperUnits',Paper_Units);
            set(gcf,'Units',Paper_Units);
            set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
            set(gcf,'Position',[0 0 PDF_Paper_Width Paper_Height]);
            display(['dataplot ',num2str(i),'...'])
            print(gcf,Image_File_Type,[Manuals_Dir,Plot_Filename])
        catch
            display(['Error: Problem with dataplot row ', num2str(i), ' (', Dataname,...
                '); check syntax of plot/save settings. Skipping case.'])
            continue
        end

    end
    clear S1 S2 style H M X Y P parameters
end

clear A

% Pack data for use in scatplot
try
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
                  {Save_Metric_Type'},...
                  {Save_Measured_Quantity},...
                  {Save_Predicted_Quantity}];
catch
    saved_data = [];
    display(['Error: saved_data = []'])
end

display('dataplot completed successfully!')


