% R. McDermott
% 7-06-2009
% scatplot.m
%
% Generate scatter plots.  Must first run dataplot.m to generate
% saved_data and drange.
%
% [] = scatplot(saved_data,drange,qfil,plotfil)
%
% Arguments:
%    saved_data - cell array of packed data, may be obtained by running
%    the dataplot.m function
%
%    drange - obtained from (or input to) function dataplot.m
%
%    qfil - file containing the scatterplot parameters
%
%    plotdir - directory where the output files are to go
%
%
% Dependencies:
%    ../scripts/define_qrow_variables.m
%

function [] = scatplot(saved_data,drange,qfil,plotdir,varargin)

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
Save_Metric_Type      = saved_data{:,11};

% If a statistics output file is specified, then enable statistics throughout
if length(varargin) >= 1
    output_file = varargin{1};
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
    output_stats{1,4} = 'Type of Metric';
    output_stats{1,5} = 'Expected Metric';
    output_stats{1,6} = 'Predicted Metric';
    output_stats{1,7} = 'Dependent Variable';
    output_stats{1,8} = 'Type of Error';
    output_stats{1,9} = 'Error';
    output_stats{1,10} = 'Error Tolerance';
    output_stats{1,11} = 'Within Specified Error Tolerance';
    output_stats{1,12} = 'Plot Filename';
    stat_line = 2;
end

for j=2:length(Q);
    
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
            size_measured = size(nonzeros(Measured_Metric(k,:,:)));
            size_predicted = size(nonzeros(Predicted_Metric(k,:,:)));
            % Skip case if predicted metric is zero
            if size_predicted(1)==0
                display(['Error: Size of predicted metric is zero for scatterplot ', Scatter_Plot_Title, '. Skipping scatterplot.'])
                continue
            end
            % Check to see if measured and predicted arrays are the same size
            if size_measured(1) ~= size_predicted(1)
                display(['Error: Mismatched measured and predicted arrays in scatter plot for scatterplot ', Scatter_Plot_Title, '. Verify that the statistical metrics are being used properly for all cases. Skipping scatterplot.'])
                continue
            end
            
            if strcmp(Plot_Type,'linear')
                K(k) = plot(nonzeros(Measured_Metric(k,:,:)),nonzeros(Predicted_Metric(k,:,:)),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
            elseif strcmp(Plot_Type,'loglog')
                K(k) = loglog(nonzeros(Measured_Metric(k,:,:)),nonzeros(Predicted_Metric(k,:,:)),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
            elseif strcmp(Plot_Type,'semilogx')
                K(k) = semilogx(nonzeros(Measured_Metric(k,:,:)),nonzeros(Predicted_Metric(k,:,:)),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
            elseif strcmp(Plot_Type,'semilogy')
                K(k) = semilogy(nonzeros(Measured_Metric(k,:,:)),nonzeros(Predicted_Metric(k,:,:)),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
            end
            
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
                        within_tolerance = 'No';
                    end
                    
                    % Write descriptive statistics to output_stats cell
                    output_stats{stat_line,1} = i;
                    output_stats{stat_line,2} = Save_Group_Key_Label{i,1};
                    output_stats{stat_line,3} = Save_Dataname{i,1};
                    output_stats{stat_line,4} = Save_Metric_Type{i,1};
                    output_stats{stat_line,5} = single_measured_metric(m);
                    output_stats{stat_line,6} = single_predicted_metric(m);
                    output_stats{stat_line,7} = Save_Dep_Title{i,1};
                    output_stats{stat_line,8} = error_type;
                    output_stats{stat_line,9} = sprintf('%1.2e', error_val);
                    output_stats{stat_line,10} = sprintf('%1.2e', error_tolerance);
                    output_stats{stat_line,11} = within_tolerance;
                    output_stats{stat_line,12} = Save_Plot_Filename{i,1};
                    stat_line = stat_line + 1;
                end
            end
        end
    end
    
    if k>0
        
        % statistics
        
        Measured_Values  = nonzeros(Measured_Metric);
        Predicted_Values = nonzeros(Predicted_Metric);
        n_pts = length(Measured_Values);
        for ib=1:10
            bin_indices = find(Measured_Values>=(ib-1)*Plot_Max/10 & Measured_Values<ib*Plot_Max/10);
            bin_weight(ib) = n_pts/length(bin_indices);
            clear bin_indices
        end
        weight = zeros(size(Measured_Values));
        for iv=1:n_pts
            for ib=1:10
                if Measured_Values(iv)>=(ib-1)*Plot_Max/10 & Measured_Values(iv)<ib*Plot_Max/10; weight(iv) = bin_weight(ib); end
            end
        end
            
        E_bar = sum(log(Measured_Values).*weight)/sum(weight);
        M_bar = sum(log(Predicted_Values).*weight)/sum(weight);
        size_measured = size(Measured_Values);
        size_predicted = size(Predicted_Values);
        % Check to see if measured and predicted arrays are the same size
        if size_measured(1) ~= size_predicted(1)
            display(['Error: Mismatched measured and predicted arrays for scatterplot ', Scatter_Plot_Title, '. Skipping scatterplot.'])
            continue
        end
        
        u2 = sum(    (((log(Predicted_Values)-log(Measured_Values)) - (M_bar-E_bar)).^2).*weight   )/(sum(weight)-1);
        u  = sqrt(u2);
        Sigma_E = Sigma_2_E/200;
        Sigma_E = min(u/sqrt(2),Sigma_E);
        Sigma_M = sqrt( max(0,u*u - Sigma_E.^2) );
        delta = exp(M_bar-E_bar+0.5*Sigma_M.^2-0.5*Sigma_E.^2);
        
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')                    
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1+2*Sigma_E)],'k--') 
        plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1-2*Sigma_E)],'k--') 
       
         if strcmp(Model_Error,'yes') 
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max],'r-')
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1+2*Sigma_M)],'r--')
             plot([Plot_Min,Plot_Max],[Plot_Min,delta*Plot_Max*(1-2*Sigma_M)],'r--')
         end
        
        % format the legend and axis labels
        xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
        ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size,'FontName',Font_Name)
        axis([Plot_Min Plot_Max Plot_Min Plot_Max])
        
        set(gca,'Units','inches')
        set(gca,'FontName','Times')
        set(gca,'FontSize',12)
        set(gca,'YTick',get(gca,'XTick'))
        set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
        
        if strcmp(Plot_Type,'linear')
            text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
            Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName','Times','Interpreter',Font_Interpreter)
        elseif strcmp(Plot_Type,'loglog')
            text(10^(log10(Plot_Min)+Title_Position(1)*(log10(Plot_Max)-log10(Plot_Min))),10^(log10(Plot_Min)+Title_Position(2)*(log10(Plot_Max)-log10(Plot_Min))),...
            Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName','Times','Interpreter',Font_Interpreter)
        elseif strcmp(Plot_Type,'semilogx')
            text(10^(log10(Plot_Min)+Title_Position(1)*(log10(Plot_Max)-log10(Plot_Min))),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
            Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName','Times','Interpreter',Font_Interpreter)
        elseif strcmp(Plot_Type,'semilogy')
            text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),10^(log10(Plot_Min)+Title_Position(2)*(log10(Plot_Max)-log10(Plot_Min))),...
            Scatter_Plot_Title,'FontSize',Scat_Title_Font_Size,'FontName','Times','Interpreter',Font_Interpreter)
        end
  
        if Sigma_E > 0.0
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.05)*(Plot_Max-Plot_Min),...
                 ['$2 \, \tilde{\sigma}_E$=',num2str(2*Sigma_E,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
        end
         
        if strcmp(Model_Error,'yes')
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.10)*(Plot_Max-Plot_Min),...
                ['$2 \, \tilde{\sigma}_M$=',num2str(2*Sigma_M,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
            
            text(Plot_Min+(Title_Position(1)+0.05)*(Plot_Max-Plot_Min),Plot_Min+(Title_Position(2)-0.15)*(Plot_Max-Plot_Min),...
                ['Bias =',num2str(delta,'%4.2f')],'FontSize',12,'FontName','Times','Interpreter',Font_Interpreter)
        end
        
        C = stripcell(Group_Key_Label);
        [B I] = unique(C);
        
        if size(Key_Position)>0
            legend_handle = legend(K(I),C(I),'Location',Key_Position,'FontSize',12','Interpreter',Font_Interpreter);
            if isequal(Key_Position,'EastOutside')
               pos = get(legend_handle,'position');
               set(legend_handle,'position',[Scat_Paper_Width pos(2:4)])
            end
            if isequal(Key_Position,'SouthEastOutside')
               pos = get(legend_handle,'position');
               set(legend_handle,'position',[Scat_Paper_Width 0.5 pos(3:4)])
            end
            set(legend_handle,'Interpreter',Font_Interpreter);
            set(legend_handle,'Fontsize',Key_Font_Size);
            set(legend_handle,'Box','on');
        end
        
        hold off
        
        % print to pdf
        
        PDF_Paper_Width = Paper_Width_Factor * Scat_Paper_Width;

        set(gcf,'Visible','on');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[PDF_Paper_Width Scat_Paper_Height]);
        set(gcf,'PaperPosition',[0 0 PDF_Paper_Width Scat_Paper_Height]);
        display(['Printing scatter plot ',num2str(j),'...'])
        print(gcf,'-dpdf',[plotdir,Plot_Filename])
        
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
        expected_value = m{i_row, 5};
        predicted_value = m{i_row, 6};
        dependent_variable = m{i_row, 7};
        % Remove " Error" from string to save horizontal space
        error_type = strrep(m{i_row, 8}, ' Error', '');
        % Convert strings to numbers for later formatting
        error_val = str2num(m{i_row, 9});
        tol = str2num(m{i_row, 10});
        % Additional columns
        within_tolerance = m{i_row, 11};
        
        % Write out all columns to .tex file
        fprintf(fid, '%s',   case_name, ' & ');
        fprintf(fid, '%s',   num2str(expected_value, '%1.2e'), ' & ');
        fprintf(fid, '%s',   num2str(predicted_value, '%1.2e'), ' & ');
        fprintf(fid, '%s',   dependent_variable, ' & ');
        fprintf(fid, '%s',   error_type, ' & ');
        fprintf(fid, '%s',   num2str(error_val, '%1.2e'), ' & ');
        fprintf(fid, '%s',   num2str(tol, '%1.2e'), ' & ');
        fprintf(fid, '%s%s\n', within_tolerance, ' \\');
    end
    fprintf(fid,'%s\n','\end{longtable}');
    fprintf(fid,'%s\n','\end{center}');
end


display('scatplot completed successfully!')


