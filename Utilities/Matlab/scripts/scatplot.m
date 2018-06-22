% R. McDermott
% 7-06-2009
% scatplot.m
%
% Generate scatter plots.  Must first run dataplot.m to generate
% saved_data and drange.
%
% [] = scatplot(saved_data, drange, Scatterplot_Inputs_File, Manuals_Dir)
%
% Arguments:
%    saved_data - cell array of packed data, may be obtained by running
%    the dataplot.m function
%
%    drange - obtained from (or input to) function dataplot.m
%
%    Scatterplot_Inputs_File - file containing the scatterplot parameters
%
%    Manuals_Dir - directory where the output files are to go
%
%
% Dependencies:
%    ../scripts/define_qrow_variables.m
%

function [] = scatplot(saved_data,drange,varargin)

% Parse input arguments using ('Key', Value) syntax
for k=1:2:length(varargin);
    switch (varargin{k})
    case {'Manuals_Dir'}
        Manuals_Dir = varargin{k+1};
    case {'Scatterplot_Inputs_File'}
        Scatterplot_Inputs_File = varargin{k+1};
    case {'Stats_Output'}
        Stats_Output = varargin{k+1};
    case {'Scatterplot_Dir'}
        Scatterplot_Dir = varargin{k+1};
    end
end

if strcmp(Stats_Output, 'Verification')
   Statistics_Tex_Output = [Scatterplot_Dir,'verification_statistics.tex'];
   Output_File = [Scatterplot_Dir,'verification_scatterplot_output.csv'];
end

if strcmp(Stats_Output, 'Validation')
   Statistics_Tex_Output = [Scatterplot_Dir,'validation_statistics.tex'];
   Output_File = [Scatterplot_Dir,'validation_scatterplot_output.csv'];
   Histogram_Tex_Output = [Scatterplot_Dir,'validation_histograms.tex'];
end

% Unpack data
Save_Quantity         = saved_data{:,1};
Save_Group_Style      = saved_data{:,2};
Save_Fill_Color       = saved_data{:,3};
Save_Group_Key_Label  = saved_data{:,4};
Save_Measured_Metric  = saved_data{:,5};
Save_Predicted_Metric = saved_data{:,6};
Save_Dataname         = saved_data{:,7};
Save_Plot_Filename    = saved_data{:,8};
Save_Dep_Title        = saved_data{:,9};
Save_Error_Tolerance  = saved_data{:,10};
Save_Metric_Type      = saved_data{:,11};
Save_Measured_Quantity= saved_data{:,12};
Save_Predicted_Quantity= saved_data{:,13};

Size_Save_Quantity = size(Save_Quantity);

% If a statistics output file is specified, then enable statistics throughout.
% This is also used to enable histogram plotting for validation cases.
% Stats_Output = 'None': No output statistics
% Stats_Output = 'Verification': FDS verification statistics
% Stats_Output = 'Validation': FDS or Correlation validation statistics
if exist('Stats_Output', 'var') == 0
    Stats_Output = 'None';
end

% Read in global plot options
plot_style

% Read in scatter plot inputs file
Q = importdata(Scatterplot_Inputs_File);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

% Generate header information for verification output_stats
if strcmp(Stats_Output, 'Verification')
    output_stats = {};
    output_stats{1,1} = 'Dataplot Line Number';
    output_stats{1,2} = 'Verification Group';
    output_stats{1,3} = 'Case Name';
    output_stats{1,4} = 'Type of Metric';
    output_stats{1,5} = 'Expected Quantity';
    output_stats{1,6} = 'Expected Value';
    output_stats{1,7} = 'Predicted Quantity';
    output_stats{1,8} = 'Predicted Value';
    output_stats{1,9} = 'Dependent Variable';
    output_stats{1,10} = 'Type of Error';
    output_stats{1,11} = 'Error';
    output_stats{1,12} = 'Error Tolerance';
    output_stats{1,13} = 'Within Specified Error Tolerance';
    output_stats{1,14} = 'Plot Filename';
    stat_line = 2;
end

% Generate header information for validation output_stats
if strcmp(Stats_Output, 'Validation')
    output_stats = {};
    output_stats{1,1} = 'Quantity';
    output_stats{1,2} = 'Number of Datasets';
    output_stats{1,3} = 'Number of Points';
    output_stats{1,4} = 'Sigma_Experiment';
    output_stats{1,5} = 'Sigma_Model';
    output_stats{1,6} = 'Bias';
    stat_line = 2;
    Output_Histograms = {};
end

file_line = strings(200);

for j=2:length(Q);

    define_qrow_variables

    Model_Error = 'yes';
    if Sigma_E < 0 ; Model_Error = 'no'; end

    clear Measured_Metric
    clear Predicted_Metric

    figure
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    set(gca,'Units','inches')
    set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

    n_lines = 0;
    k = 0;
    for i=drange
        if i > Size_Save_Quantity(2); break; end
        if strcmp(Save_Quantity(1,i),Scatter_Plot_Title) || strcmp(Save_Quantity(Size_Save_Quantity(1),i),Scatter_Plot_Title)
            k = k+1;

            Measured_Metric(k,:,:)  = Save_Measured_Metric(i,:,:);
            Predicted_Metric(k,:,:) = Save_Predicted_Metric(i,:,:);

            for kk=1:10
                for jj=1:10
                    if abs(Measured_Metric(k,jj,kk))>0
                       n_lines = n_lines+1;
                       file_line(n_lines) = [num2str(i),',',num2str(Measured_Metric(k,jj,kk)),',',num2str(Predicted_Metric(k,jj,kk))];
                    end
                    if Measured_Metric(k,jj,kk) < Plot_Min; Measured_Metric(k,jj,kk)=0. ; Predicted_Metric(k,jj,kk)=0. ; end
                    if Measured_Metric(k,jj,kk) > Plot_Max; Measured_Metric(k,jj,kk)=0. ; Predicted_Metric(k,jj,kk)=0. ; end
                    if Predicted_Metric(k,jj,kk) < Plot_Min; Measured_Metric(k,jj,kk)=0. ; Predicted_Metric(k,jj,kk)=0. ; end
                    if Predicted_Metric(k,jj,kk) > Plot_Max; Measured_Metric(k,jj,kk)=0. ; Predicted_Metric(k,jj,kk)=0. ; end
                end
            end

            Nonzeros_Measured_Metric = nonzeros(Measured_Metric(k,:,:));
            Nonzeros_Predicted_Metric = nonzeros(Predicted_Metric(k,:,:));
            Size_Measured = size(Nonzeros_Measured_Metric);
            Size_Predicted = size(Nonzeros_Predicted_Metric);
            if Size_Measured(1) == 0 ; k=k-1 ; continue ; end
            Group_Key_Label(k) = Save_Group_Key_Label(i);

            % Check to see if measured and predicted arrays are the same size
            if Size_Measured(1) ~= Size_Predicted(1)
                display(['Error: Mismatched measured and predicted arrays in scatter plot for scatterplot ', Scatter_Plot_Title, '. Verify that the statistical metrics are being used properly for all cases. Skipping scatterplot.'])
                continue
            end

            if strcmp(Plot_Type,'linear')
                K(k) = plot(Nonzeros_Measured_Metric,Nonzeros_Predicted_Metric,...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i)),'MarkerSize',Marker_Size); hold on
            elseif strcmp(Plot_Type,'loglog')
                K(k) = loglog(Nonzeros_Measured_Metric,Nonzeros_Predicted_Metric,...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i)),'MarkerSize',Marker_Size); hold on
            elseif strcmp(Plot_Type,'semilogx')
                K(k) = semilogx(Nonzeros_Measured_Metric,Nonzeros_Predicted_Metric,...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i)),'MarkerSize',Marker_Size); hold on
            elseif strcmp(Plot_Type,'semilogy')
                K(k) = semilogy(Nonzeros_Measured_Metric,Nonzeros_Predicted_Metric,...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i)),'MarkerSize',Marker_Size); hold on
            end

            % Perform this code block for FDS verification scatterplot output
            if strcmp(Stats_Output, 'Verification')
                single_measured_metric = Nonzeros_Measured_Metric;
                single_predicted_metric = Nonzeros_Predicted_Metric;
                % Loop over multiple line comparisons and build output_stats cell
                for m=1:length(single_measured_metric)

                    % Get type of statistics to compute
                    error_type = Save_Quantity{1,i};

                    % Compute the appropriate type of statistics, depending
                    % on the 'Quantity' specification in dataplot_inputs
                    if strcmp(error_type, 'Relative Error')
                        error_val = abs((single_predicted_metric(m)-single_measured_metric(m))/single_measured_metric(m));
                    elseif strcmp(error_type, 'Absolute Error')
                        error_val = abs(single_predicted_metric(m)-single_measured_metric(m));
                    end

                    % Compare the error to the specified error tolerance,
                    % which is in the dataplot_inputs column called 'Error_Tolerance'
                    error_tolerance = str2num(Save_Error_Tolerance{i,1});
                    if error_val <= error_tolerance
                        within_tolerance = 'Yes';
                    else
                        within_tolerance = 'Out of Tolerance';
                    end

                    % Write verification statistics to output_stats cell
                    output_stats{stat_line,1} = i;
                    output_stats{stat_line,2} = Save_Group_Key_Label{i,1};
                    output_stats{stat_line,3} = Save_Dataname{i,1};
                    output_stats{stat_line,4} = Save_Metric_Type{i,1};
                    output_stats{stat_line,5} = Save_Measured_Quantity{i,m,1};
                    output_stats{stat_line,6} = single_measured_metric(m);
                    output_stats{stat_line,7} = Save_Predicted_Quantity{i,m,1};
                    output_stats{stat_line,8} = single_predicted_metric(m);
                    output_stats{stat_line,9} = Save_Dep_Title{i,1};
                    output_stats{stat_line,10} = error_type;
                    output_stats{stat_line,11} = sprintf('%1.2e', error_val);
                    output_stats{stat_line,12} = sprintf('%1.2e', error_tolerance);
                    output_stats{stat_line,13} = within_tolerance;
                    output_stats{stat_line,14} = Save_Plot_Filename{i,1};
                    stat_line = stat_line + 1;
                end
            end
        end
    end

    % Write out the measured and predicted value for this quantity
    if n_lines > 0 && strcmp(Stats_Output, 'Validation')
       fid = fopen([Manuals_Dir,Plot_Filename,'_Points.csv'],'w');
       fprintf(fid,'%s\n','Row,Measured Value,Predicted Value');
       for nn=1:n_lines
          fprintf(fid,'%s\n',file_line(nn));
       end
       fclose(fid);
    end

    % Perform this code block for FDS validation scatterplot output
    if k > 0 && strcmp(Stats_Output, 'Validation')

        Measured_Values = nonzeros(Measured_Metric);
        Predicted_Values = nonzeros(Predicted_Metric);
        n_pts = length(Measured_Values);

        % Weight the data -- for each point on the scatterplot compute a
        % "weight" to provide sparse data with greater importance in the
        % calculation of the accuracy statistics

        weight = zeros(size(Measured_Values));

        if strcmp(Weight_Data,'yes')
            Max_Measured_Value = max(Measured_Values);
            Bin_Size = Max_Measured_Value/10;
            for ib=1:10
                bin_indices = find(Measured_Values>(ib-1)*Bin_Size & Measured_Values<=ib*Bin_Size);
                bin_weight(ib) = n_pts/length(bin_indices);
                clear bin_indices
            end
            for iv=1:n_pts
                for ib=1:10
                    if Measured_Values(iv)>(ib-1)*Bin_Size && Measured_Values(iv)<=ib*Bin_Size; weight(iv) = bin_weight(ib); end
                end
            end
        else
            for iv=1:n_pts
                weight(iv) = 1.;
            end
        end

        % Calculate statistics

        log_E_bar = sum(log(Measured_Values).*weight)/sum(weight);
        log_M_bar = sum(log(Predicted_Values).*weight)/sum(weight);
        Size_Measured = size(Measured_Values);
        Size_Predicted = size(Predicted_Values);

        if Size_Measured(1) ~= Size_Predicted(1)
            display(['Error: Mismatched measured and predicted arrays for scatterplot ', Scatter_Plot_Title, '. Skipping scatterplot.'])
            continue
        end

        u2 = sum(    (((log(Predicted_Values)-log(Measured_Values)) - (log_M_bar-log_E_bar)).^2).*weight   )/(sum(weight)-1);
        u  = sqrt(u2);
        Sigma_E = Sigma_E/100;
        Sigma_E = min(u/sqrt(2),Sigma_E);
        Sigma_M = sqrt( max(0,u*u - Sigma_E.^2) );
        delta = exp(log_M_bar-log_E_bar+0.5*Sigma_M.^2-0.5*Sigma_E.^2);

        % Plot diagonal lines
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')
        if strcmp(Model_Error, 'yes')
            plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max]*(1+2*Sigma_E),'k--')
            plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max]/(1+2*Sigma_E),'k--')
            plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max]*delta,'r-')
         %  if delta > 2*Sigma_M
               plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max]*delta*(1+2*Sigma_M),'r--')
               plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max]*delta/(1+2*Sigma_M),'r--')
         %  end
        end

        % Format the legend and axis labels
        xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
        ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
        axis([Plot_Min Plot_Max Plot_Min Plot_Max])

        set(gca,'FontName',Font_Name)
        set(gca,'FontSize',Key_Font_Size)
        set(gca,'YTick',get(gca,'XTick'))

        curtick = get(gca, 'XTick');
        set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
        curtick = get(gca, 'YTick');
        set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));

        if strcmp(Plot_Type,'linear') | strcmp(Plot_Type,'semilogy')
            Title_Position_X = Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min);
        else
            Title_Position_X = 10^(log10(Plot_Min)+Title_Position(1)*(log10(Plot_Max)-log10(Plot_Min)));
        end

        if strcmp(Plot_Type,'linear') | strcmp(Plot_Type,'semilogx')
            Title_Position_Y_1 = Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min);
            Title_Position_Y_2 = Plot_Min+(Title_Position(2)-0.05)*(Plot_Max-Plot_Min);
            Title_Position_Y_3 = Plot_Min+(Title_Position(2)-0.10)*(Plot_Max-Plot_Min);
            Title_Position_Y_4 = Plot_Min+(Title_Position(2)-0.15)*(Plot_Max-Plot_Min);
        else
            Title_Position_Y_1 = 10^(log10(Plot_Min)+Title_Position(2)*(log10(Plot_Max)-log10(Plot_Min)));
            Title_Position_Y_2 = 10^(log10(Plot_Min)+(Title_Position(2)-0.05)*(log10(Plot_Max)-log10(Plot_Min)));
            Title_Position_Y_3 = 10^(log10(Plot_Min)+(Title_Position(2)-0.10)*(log10(Plot_Max)-log10(Plot_Min)));
            Title_Position_Y_4 = 10^(log10(Plot_Min)+(Title_Position(2)-0.15)*(log10(Plot_Max)-log10(Plot_Min)));
        end

        text(Title_Position_X,Title_Position_Y_1,Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)

        if Sigma_E > 0.0
            text(Title_Position_X,Title_Position_Y_2,['Exp. Rel. Std. Dev.: ',num2str(Sigma_E,'%4.2f')],'FontSize',Key_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
        end

        if strcmp(Model_Error,'yes')
            text(Title_Position_X,Title_Position_Y_3,['Model Rel. Std. Dev.: ',num2str(Sigma_M,'%4.2f')],'FontSize',Key_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
            text(Title_Position_X,Title_Position_Y_4,['Model Bias Factor: ',num2str(delta,'%4.2f')],'FontSize',Key_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
        end

        C = stripcell(Group_Key_Label);
        [B I] = unique(C);

        if size(Key_Position) > 0
            legend_handle = legend(K(I),C(I),'Location',Key_Position,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter,'Units',Plot_Units);
            if strcmp(Key_Position,'EastOutside')
               pos = get(legend_handle,'position');
               set(legend_handle,'position',[Scat_Paper_Width pos(2:4)])
            end
            if strcmp(Key_Position,'SouthEastOutside')
               pos = get(legend_handle,'position');
               set(legend_handle,'position',[Scat_Paper_Width pos(2:4)])
            end
            set(legend_handle,'Interpreter',Font_Interpreter);
            set(legend_handle,'Fontsize',Key_Font_Size);
            set(legend_handle,'Box','on');
        end

        hold off

        % Print to pdf
        PDF_Paper_Width = Paper_Width_Factor * Scat_Paper_Width;

        set(gcf,'Visible',Figure_Visibility);
        set(gcf,'Units',Paper_Units);
        set(gcf,'PaperUnits',Paper_Units);
        set(gcf,'PaperSize',[PDF_Paper_Width Scat_Paper_Height]);
        set(gcf,'Position',[0 0 PDF_Paper_Width Scat_Paper_Height]);
        display(['scatplot ',num2str(j),'...'])
        print(gcf,Image_File_Type,[Manuals_Dir,Plot_Filename])

        % Print histogram of ln(M/E) and normal distribution
        statistics_histogram

        % Write validation statistics to output_stats cell
        output_stats{stat_line,1} = Scatter_Plot_Title; % Quantity
        output_stats{stat_line,2} = size(B, 2); % Number of data sets
        output_stats{stat_line,3} = size(Predicted_Values, 1); % Number of data points
        output_stats{stat_line,4} = sprintf('%0.2f', Sigma_E); % Sigma_E
        output_stats{stat_line,5} = sprintf('%0.2f', Sigma_M); % Sigma_M
        output_stats{stat_line,6} = sprintf('%0.2f', delta); % Bias
        stat_line = stat_line + 1;

    elseif k == 0 && strcmp(Stats_Output, 'Validation')
        display(['No data for scatter plot ', Scatter_Plot_Title])
    end

    clear Measured_Metric Measured_Values Predicted_Metric Predicted_Values Group_Key_Label K
    close all
end

% Verification and validation statistics output routines
statistics_output

fclose('all');

display('scatplot completed successfully!')


