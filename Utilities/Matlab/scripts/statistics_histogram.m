% K. Overholt
% 8-11-2013
% statistics_histogram.m

% This script is called from scatplot and writes out
% histogram plots to be used in various guides.

% Print histogram of ln(M/E) and normal distribution
if strcmp(Stats_Output, 'Validation')
    % Wrap histogram routine in try loop
    % Skips case upon any Matlab error
    try
        ln_M_E = log(Predicted_Values)-log(Measured_Values);
        % Normality test (requires at least 4 observations)
        if length(ln_M_E) >= 4

            pval = spiegel_test(ln_M_E);

            % Plot histogram
            figure
            box on
            hold on
            [n,xout] = hist(ln_M_E,10);
            bar(xout,n,1,'LineWidth',1,'FaceColor',[0.7,0.7,0.7])

            % Plot normal distribution
            x_lim = [xout(1)-(xout(2)-xout(1)),xout(end)+(xout(2)-xout(1))];
            ix = x_lim(1):1e-3:x_lim(2);
            mu = mean(ln_M_E);
            sd = std(ln_M_E);
            % Generate PDF of normal distribution (overlaid on results)
            iy = 1/(sd*sqrt(2*pi))*exp(-(ix-mu).^2/(2*sd^2));
            plot(ix,iy*trapz(xout,n),'k','LineWidth',2);

            % Additional plot content
            set(gca,'XLim',[x_lim(1),x_lim(2)]);
            y_lim = get(gca,'YLim') * 1.25;
            set(gca,'YLim',y_lim);
            xlabel('Interval Number','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
            ylabel('Number of Data Points','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
            set(gca,'FontName',Font_Name)
            set(gca,'FontSize',12)
            set(gca,'XTick',xout,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'})
            text(0.03, 0.90,Scatter_Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter,'Units','normalized')
            text(0.03, 0.82,'Normality Test','FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter,'Units','normalized')
            text(0.03, 0.74,['p-value = ',num2str(pval,'%4.2f')],'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter,'Units','normalized')

            PDF_Paper_Width = Paper_Width;

            set(gcf,'Visible',Figure_Visibility);
            set(gcf,'Units',Paper_Units);
            set(gcf,'PaperUnits',Paper_Units);
            set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
            set(gcf,'Position',[0 0 PDF_Paper_Width Paper_Height]);
            print(gcf,Image_File_Type,[Manuals_Dir,[Plot_Filename, '_Histogram']])
            hold off
            % Add histogram name to array for LaTeX output later
            [~, filename, ~] = fileparts(Plot_Filename);
            Output_Histograms{end+1} = [filename, '_Histogram'];
        end
    catch
        display(['Error: Problem with histogram routine for scatter plot ', Scatter_Plot_Title, '; Skipping histogram.'])
    end
end
