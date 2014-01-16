% K. Overholt
% 8-7-2013
% statistics_output.m

% This script is called from scatplot and writes out statistical
% information (csv and tex files) to be used in various guides.

% Write all verification or validation statistics from output_stats to csv output_file
if ~strcmp(Stats_Output, 'None')
    [rows, cols] = size(output_stats);
    fid = fopen(Output_File, 'w');
    for i_row = 1:rows
        file_line = '';
        for i_col = 1:cols
            contents = output_stats{i_row, i_col};
            if isnumeric(contents)
                contents = num2str(contents);
            elseif isstr(contents)
                contents = strcat('"', contents, '"');
            elseif isempty(contents)
                contents = '';
            end
            if i_col < cols
                file_line = [file_line, contents, ','];
            else
                file_line = [file_line, contents];
            end
        end
        count = fprintf(fid, '%s\n', file_line);
    end
    fclose(fid);
end

% Write statistics information to a LaTeX table
% for inclusion in the FDS Verification Guide
if strcmp(Stats_Output, 'Verification')
    fid = fopen(Statistics_Tex_Output, 'wt');
    % Generate table header information in .tex file
    fprintf(fid, '%s\n', '\tiny');
    fprintf(fid, '%s\n', '\begin{longtable}[c]{|l|c|c|c|c|c|c|}');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', 'Case Name & Expected & Predicted & Type of Error & Error & Error     & Within    \\');
    fprintf(fid, '%s\n', '          & Metric   & Metric    &               &       & Tolerance & Tolerance \\ \hline \hline');
    fprintf(fid, '%s\n', '\endfirsthead');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', 'Case Name & Expected & Predicted & Type of Error & Error & Error     & Within    \\');
    fprintf(fid, '%s\n', '          & Metric   & Metric    &               &       & Tolerance & Tolerance \\ \hline \hline');
    fprintf(fid, '%s\n', '\endhead');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', '\endfoot');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', '\endlastfoot');
    [rows, cols] = size(output_stats);
    for i_row = 2:rows
        % Format strings for various columns in table (and add short names)
        m = output_stats;
        % Escape underscores for LaTeX
        case_name = strrep(m{i_row, 3}, '_', '\_');
        % Additional columns
        expected_value = m{i_row, 5};
        predicted_value = m{i_row, 6};
        % Remove " Error" from string to save horizontal space
        error_type = strrep(m{i_row, 8}, ' Error', '');
        % Convert strings to numbers for later formatting
        error_val = str2num(m{i_row, 9});
        tol = str2num(m{i_row, 10});
        % Additional columns
        within_tolerance = m{i_row, 11};
        
        % Write out all columns to .tex file
        fprintf(fid, '%s', case_name, ' & ');
        fprintf(fid, '%s', num2str(expected_value, '%1.2e'), ' & ');
        fprintf(fid, '%s', num2str(predicted_value, '%1.2e'), ' & ');
        fprintf(fid, '%s', error_type, ' & ');
        fprintf(fid, '%s', num2str(error_val, '%1.2e'), ' & ');
        fprintf(fid, '%s', num2str(tol, '%1.2e'), ' & ');
        fprintf(fid, '%s%s\n', within_tolerance, ' \\');
    end
    fprintf(fid,'%s\n','\end{longtable}');
    fprintf(fid, '%s\n', '\normalsize');
end

% Write statistics information to a LaTeX table for inclusion
% in the FDS Validation Guide or Correlation Guide
if strcmp(Stats_Output, 'Validation')
    fid = fopen(Statistics_Tex_Output, 'wt');
    % Generate table header information in .tex file
    fprintf(fid, '%s\n', '\begin{longtable}[c]{|l|c|c|c|c|c|}');
    fprintf(fid, '%s\n', '\caption[Summary statistics]{Summary statistics for all quantities of interest}');
    fprintf(fid, '%s\n', '\label{summary_stats}');
    fprintf(fid, '%s\n', '\\ \hline');
    fprintf(fid, '%s\n', 'Quantity & Datasets  & Points    & $\widetilde{\sigma}_{\rm E}$ & $\widetilde{\sigma}_{\rm M}$ & Bias \\ \hline \hline');
    fprintf(fid, '%s\n', '\endfirsthead');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', 'Quantity & Datasets  & Points    & $\widetilde{\sigma}_{\rm E}$ & $\widetilde{\sigma}_{\rm M}$ & Bias \\ \hline \hline');
    fprintf(fid, '%s\n', '\endhead');
    [rows, cols] = size(output_stats);
    for i_row = 2:rows
        % Format strings for various columns in table (and add short names)
        m = output_stats;
        quantity = m{i_row, 1};
        number_datasets = m{i_row, 2};
        number_points= m{i_row, 3};
        sigma_e = m{i_row, 4};
        sigma_m = m{i_row, 5};
        bias = m{i_row, 6};
        
        % Do not print rows with no exp. error (specified as "-1" in scatplot_inputs)
        if str2num(sigma_e) < 0
            continue
        end
        
        % Write out all columns to .tex file
        fprintf(fid, '%s', quantity, ' & ');
        fprintf(fid, '%s', num2str(number_datasets), ' & ');
        fprintf(fid, '%s', num2str(number_points), ' & ');
        fprintf(fid, '%s', num2str(sigma_e, '%0.2f'), ' & ');
        fprintf(fid, '%s', num2str(sigma_m, '%0.2f'), ' & ');
        fprintf(fid, '%s%s\n', num2str(bias, '%0.2f'), ' \\ \hline');
    end
    fprintf(fid,'%s\n','\end{longtable}');
end

% Write histogram information to a LaTeX file for inclusion
% in the FDS Validation Guide or Correlation Guide
if strcmp(Stats_Output, 'Validation') && (exist('Output_Histograms','var') == 1) && (isempty(Output_Histograms) == 0)
    fid = fopen(Histogram_Tex_Output, 'wt');
    % Write plots to LaTeX figures, eight per page
    num_histograms = length(Output_Histograms);
    page_count = ceil(num_histograms/8);
    for i = 1:page_count
        fprintf(fid, '%s\n', '\begin{figure}[p]');
        fprintf(fid, '%s\n', '\begin{tabular*}{\textwidth}{l@{\extracolsep{\fill}}r}');
        % Indices go from 1:8, 9:16, 17:24, up to last page,
        % which might have less than 8 plots
        for j = ((i-1)*8+1):(min((i*8),num_histograms-((i-1)*8)))
            % Alternate line endings in LaTeX plots
            if mod(j,2) == 1
                line_ending = '&';
            else
                line_ending = '\\';
            end
            fprintf(fid, '%s\n', ['\includegraphics[height=2.2in]{SCRIPT_FIGURES/ScatterPlots/',Output_Histograms{j},'} ',line_ending]);
        end
        fprintf(fid, '%s\n', '\end{tabular*}');
        fprintf(fid, '%s\n', ['\label{Histogram_',num2str(i),'}']);
        fprintf(fid, '%s\n\n', '\end{figure}');
    end
end
