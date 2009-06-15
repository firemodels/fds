% R. McDermott
% 6-11-2009
% scatterplot.m
%
% Generate scatter plots.  Must first run read_dline.m to generate
% measured_metric and predicted_metric.
%
% Dependencies:
%    define_qrow_variables.m
%    ../scatterplot_config_matlab.csv

close all

qfil = ['../scatterplot_config_matlab.csv'];
qrange = 2:2;

addpath('../functions')
paper_width  = 6.0; % inches
paper_height = 4.0; % inches

Q = importdata(qfil);
H = textscan(Q{1},'%q','delimiter',',');
headers = H{:}'; clear H

for j=qrange
    if j>length(Q); break; end
    
    define_qrow_variables
    
    k = 1;
    for i=drange
        if strcmp(Save_Quantity(i),Scatter_Plot_Title)
            Measured_Metric(k)  = Save_Measured_Metric(i);
            Predicted_Metric(k) = Save_Predicted_Metric(i);
            K(i) = plot(Measured_Metric(k),Predicted_Metric(k),...
                char(Save_Group_Style(i)),'MarkerFaceColor',char(Save_Fill_Color(i))); hold on
            k = k+1;
        end
    end
    
    % statistics
    plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max],'k-')                    % predicted = measured
    plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1+Sigma_2_E/100)],'k--') % + Sigma_2_E
    plot([Plot_Min,Plot_Max],[Plot_Min,Plot_Max*(1-Sigma_2_E/100)],'k--') % - Sigma_2_E
    
    % format the legend and axis labels
    xlabel(Ind_Title,'Interpreter','LaTeX','FontSize',14)
    ylabel(Dep_Title,'Interpreter','LaTeX','FontSize',14)
    axis([Plot_Min Plot_Max Plot_Min Plot_Max])
    
    set(gca,'FontName','Times')
    set(gca,'FontSize',12)
    set(gca,'YTick',get(gca,'XTick'))
    
    text(Plot_Min+Title_Position(1)*(Plot_Max-Plot_Min),Plot_Min+Title_Position(2)*(Plot_Max-Plot_Min),...
        Scatter_Plot_Title,'FontSize',14,'FontName','Times','Interpreter','LaTeX')
    
    C = stripcell(Save_Group_Key_Label);
    [B I] = unique(C);
    legend(K(I),C(I),'Location',Key_Position,'FontSize',12)
    legend boxon
    
    hold off
    clear Measured_Metric Predicted_Metric
    
    % print to pdf
    set(gcf,'Visible','on');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[paper_width paper_height]);
    set(gcf,'PaperPosition',[0 0 paper_width paper_height]); 
    display(['Printing scatter plot ',num2str(j),'...'])
    print(gcf,'-dpdf',['../../../Manuals/',Plot_Filename])
end




