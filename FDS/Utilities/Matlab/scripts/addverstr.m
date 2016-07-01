% McDermott
% 6-23-2015
% addverstr.m

function []=addverstr(ha,fn,pt,varargin)

Font_Name = 'Times';
Font_Interpreter = 'TeX';

if length(varargin)==0
    VerStr_Scale_X = 0.65;
    VerStr_Scale_Y = 1.05;
else
    VerStr_Scale_X = varargin{1};
    VerStr_Scale_Y = varargin{2};
end

if exist(fn,'file')
    VerStr = importdata(fn);
    x_lim = get(ha,'XLim');
    y_lim = get(ha,'YLim');
    if strcmp(pt,'loglog')
        X_VerStr_Position = 10^( log10(x_lim(1))+ VerStr_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
        Y_VerStr_Position = 10^( log10(y_lim(1))+ VerStr_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    elseif strcmp(pt,'semilogx')
        X_VerStr_Position = 10^( log10(x_lim(1))+ VerStr_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
        Y_VerStr_Position = y_lim(1)+VerStr_Scale_Y*(y_lim(2)-y_lim(1));
    elseif strcmp(pt,'semilogy')
        X_VerStr_Position = x_lim(1)+VerStr_Scale_X*(x_lim(2)-x_lim(1));
        Y_VerStr_Position = 10^( log10(y_lim(1))+ VerStr_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    else
        X_VerStr_Position = x_lim(1)+VerStr_Scale_X*(x_lim(2)-x_lim(1));
        Y_VerStr_Position = y_lim(1)+VerStr_Scale_Y*(y_lim(2)-y_lim(1));
    end
    if isnumeric(VerStr)
        text(X_VerStr_Position,Y_VerStr_Position,['VerStr ',num2str(VerStr)], ...
            'FontSize',10,'FontName',Font_Name,'Interpreter','none')
    elseif ischar(VerStr{1})
        text(X_VerStr_Position,Y_VerStr_Position,[VerStr], ...
            'FontSize',10,'FontName',Font_Name,'Interpreter','none')
    end
end