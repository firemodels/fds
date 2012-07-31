% R. McDermott
% 7-06-2009
% scatplot.m
%
% Generate scatter plots.  Must first run dataplot.m to generate
% saved_data and drange.
%
% [] = scatplot(saved_data,drange,[qrange])
%
% Arguments:
%    saved_data - cell array of packed data, may be obtained by running
%    the dataplot.m function
%
%    drange - obtained from (or input to) function dataplot.m
%
%    [optional] qrange - specify the range of lines to read from the
%    scatterplot configuration file
%
% Dependencies:
%    ../scripts/define_qrow_variables.m
%
% Example:
%    From the Matlab/functions/ directory, type
%    >> scatplot(saved_data,drange,3:4)

function [] = scatplot(saved_data,drange,varargin)

% unpack data
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

qfil = varargin{1};

% If a statistics output file is specified, then enable statistics throughout
if length(varargin) >= 2
    output_file = varargin{2};
    stats_output = 1;
else
    stats_output = 0;
end

qrange = [2:100];

plot_style

Q = importdata(qfil);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

% Generate header information for output_stats
if stats_output == 1
    output_stats = {};
    output_stats{1,1} = 'Dataplot Line Number';
    output_stats{1,2} = 'Verification Group';
    output_stats{1,3} = 'Case Name';
    output_stats{1,4} = 'Expected Metric';
    output_stats{1,5} = 'Predicted Metric';
    output_stats{1,6} = 'Dependent Variable';
    output_stats{1,7} = 'Type of Error';
    output_stats{1,8} = 'Error';
    output_stats{1,9} = 'Error Tolerance';
    output_stats{1,10} = 'Within Specified Error Tolerance';
    output_stats{1,11} = 'Plot Filename';
    stat_line = 2;
end

for j=qrange
    if j>length(Q); break; end
    
    define_qrow_variables
    
    clear Measured_Metric
    clear Predicted_Metric
    
    k = 0;
    for i=drange
        if i>length(Save_Quantity); break; end
        if strcmp(Save_Quantity(i),Scatter_Plot_Title)
            k = k+1;
            Measured_Metric(k,:,:)  = Save_Measured_Metric(i,:,:);
            Predicted_Metric(k,:,:) = Save_Predicted_Metric(i,:,:);
            Group_Key_Label(k)  = Save_Group_Key_Label(i);
            K(k) = plot(nonzeros(Measured_Metric(k,:,:)),nonzeros(Predicted_Metric(k,:,:)),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on

            if stats_output == 1
                single_measured_metric = nonzeros(Measured_Metric(k,:,:));
                single_predicted_metric = nonzeros(Predicted_Metric(k,:,:));
                % Loop over multiple line comparisons and build output_stats cell
                for m=1:length(single_measured_metric)
                    
                    % Get type of statistics to compute
                    error_type = Save_Quantity{i,1};
                    
                    % Compute the appropriate type of statistics, depending
                    % on the 'Quantity' specification in dataplot_inputs
                    if strcmp(error_type, 'Relative Error')
                        error = ((single_predicted_metric(m)-single_measured_metric(m))/single_measured_metric(m));
                    elseif strcmp(error_type, 'Absolute Error')
                        error = (single_predicted_metric(m)-single_measured_metric(m));
                    end
                    
                    % Compare the error to the specified error tolerance,
                    % which is in the dataplot_inputs column called 'Error_Tolerance'
                    error_tolerance = str2num(Save_Error_Tolerance{i,1});
                    if abs(error) <= error_tolerance
                        within_tolerance = 'Yes';
                    else
                        within_tolerance = 'No';
                    end
                    
                    % Write descriptive statistics to output_stats cell
                    output_stats{stat_line,1} = i;
                    output_stats{stat_line,2} = Save_Group_Key_Label{i,1};
                    output_stats{stat_line,3} = Save_Dataname{i,1};
                    output_stats{stat_line,4} = single_measured_metric(m);
                    output_stats{stat_line,5} = single_predicted_metric(m);
                    output_stats{stat_line,6} = Save_Dep_Title{i,1};
                    output_stats{stat_line,7} = error_type;
                    output_stats{stat_line,8} = sprintf('%1.8f', error);
                    output_stats{stat_line,9} = sprintf('%1.8f', error_tolerance);
                    output_stats{stat_line,10} = within_tolerance;
                    output_stats{stat_line,11} = Save_Plot_Filename{i,1};
                    stat_line = stat_line + 1;
                end
            end
        end
    end
    
    if k>0
        
        % statistics
        
        n_pts = length(nonzeros(Measured_Metric));
        E_bar = mean(log(nonzeros(Measured_Metric)));
        M_bar = mean(log(nonzeros(Predicted_Metric)));
      %  [normality,p] = lillietest(log(nonzeros(Predicted_Metric))-log(nonzeros(Measured_Metric)));
        normality = 0;
        u  = sqrt( sum( ( (log(nonzeros(Predicted_Metric))-log(nonzeros(Measured_Metric))) - (M_bar-E_bar) ).^2 )/(n_pts-1) );
        Sigma_E = Sigma_2_E/200;
        Sigma_E = min(u/sqrt(2),Sigma_E);
        Sigma_M = sqrt( max(0,u*u - Sigma_E.^2) );
        delta = exp(M_bar-E_bar+0.5*Sigma_M.^2-0.5*Sigma_E.^2);
        
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')                    
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1+2*Sigma_E)],'k--') 
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1-2*Sigma_E)],'k--') 
       
         if strcmp(Model_Error,'yes') & normality==0
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max],'r-')
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1+2*Sigma_M)],'r--')
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1-2*Sigma_M)],'r--')
         end
        
        % format the legend and axis labels
        xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
        ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
        axis([Plot_Min Plot_Max Plot_Min Plot_Max])
        
        set(gca,'Units','inches')
        set(gca,'FontName','Times')
        set(gca,'FontSize',12)
        set(gca,'YTick',get(gca,'XTick'))
        set(gca,'Position',[1,1,4.5,4.5])
        
        text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
            Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName','Times','Interpreter',Font_Interpreter)
  
         if Sigma_E > 0.0
             text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.05)*(Plot_Max-Plot_Min),...
                 ['$2 \, \tilde{\sigma}_E$=',num2str(2*Sigma_E,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
         end
         
        if strcmp(Model_Error,'yes') & normality==0
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.10)*(Plot_Max-Plot_Min),...
                ['$2 \, \tilde{\sigma}_M$=',num2str(2*Sigma_M,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
            
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.15)*(Plot_Max-Plot_Min),...
                ['Bias =',num2str(delta,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
        end
        
        C = stripcell(Group_Key_Label);
        [B I] = unique(C);
        legend(K(I),C(I),'Location',Key_Position,'FontSize',12','Interpreter',Font_Interpreter)
        legend boxon
        
        hold off
        
        % print to pdf
        set(gcf,'Visible','on');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
        set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
        display(['Printing scatter plot ',num2str(j),'...'])
        print(gcf,'-dpdf',[pwd,'/../../Manuals/',Plot_Filename])
        
    else
        display(['No data for scatter plot ',Scatter_Plot_Title])
    end
    
    clear Measured_Metric Predicted_Metric Group_Key_Label K
end

% Write all statistics from output_stats to csv output_file
if stats_output == 1
    [rows, cols] = size(output_stats);
    fid = fopen(output_file, 'w');
    for i_row = 1:rows
        file_line = '';
        for i_col = 1:cols
            contents = output_stats{i_row, i_col};
            if isnumeric(contents)
                contents = num2str(contents);
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
    st = fclose(fid);
end    

% Write statistics information to a LaTeX table for inclusion in the
% FDS Verification Guide (SCRIPT_FIGURES/verification_statistics.tex)
if stats_output == 1
    filename = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/verification_statistics.tex';
    fid = fopen(filename, 'wt');
    % Generate table header information in .tex file
    fprintf(fid,'%s\n','\begin{center}');
    fprintf(fid,'%s\n','\tiny');
    fprintf(fid,'%s\n','\begin{longtable}{|c|c|c|c|c|c|c|c|}');
    fprintf(fid,'%s\n','\hline');
    fprintf(fid, '%s\n', 'Case Name & Expected & Predicted & Dependent & Type of Error & Error & Error     & Within   \\');
    fprintf(fid, '%s\n', '          & Metric   & Metric    & Variable  &               &       & Tolerance & Tolerance \\ \hline');
    fprintf(fid, '%s\n', '\endfirsthead');
    fprintf(fid, '%s\n', '\hline');
    fprintf(fid, '%s\n', 'Case Name & Expected & Predicted & Dependent & Type of Error & Error & Error     & Within   \\');
    fprintf(fid, '%s\n', '          & Metric   & Metric    & Variable  &               &       & Tolerance & Tolerance \\ \hline');
    fprintf(fid, '%s\n', '\endhead');
    fprintf(fid, '%s\n',  '\hline');
    fprintf(fid, '%s\n', '\endfoot');
    fprintf(fid, '%s\n',  '\hline');
    fprintf(fid, '%s\n', '\endlastfoot');
    [rows, cols] = size(output_stats);
    for i_row = 2:rows
        % Format strings for various columns in table (and add short names)
        m = output_stats;
        % Escape underscores for LaTeX
        case_name = strrep(m{i_row, 3}, '_', '\_');
        % Additional columns
        expected_value = m{i_row, 4};
        predicted_value = m{i_row, 5};
        dependent_variable = m{i_row, 6};
        % Remove " Error" from string to save horizontal space
        error_type = strrep(m{i_row, 7}, ' Error', '');
        % Convert strings to numbers for later formatting
        error = str2num(m{i_row, 8});
        tol = str2num(m{i_row, 9});
        % Additional columns
        within_tolerance = m{i_row, 10};
        
        % Write out all columns to .tex file
        fprintf(fid, '%s',   case_name, ' & ');
        fprintf(fid, '%s',   num2str(expected_value, '%1.2e'), ' & ');
        fprintf(fid, '%s',   num2str(predicted_value, '%1.2e'), ' & ');
        fprintf(fid, '%s',   dependent_variable, ' & ');
        fprintf(fid, '%s',   error_type, ' & ');
        fprintf(fid, '%s',   num2str(error, '%1.2e'), ' & ');
        fprintf(fid, '%s',   num2str(tol, '%1.2e'), ' & ');
        fprintf(fid, '%s%s\n', within_tolerance, ' \\');
    end
    fprintf(fid,'%s\n','\end{longtable}');
    fprintf(fid,'%s\n','\end{center}');
end


display('scatplot completed successfully!')


