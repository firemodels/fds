% K. McGrattan
% 7-16-2009
% plot_style.m
%
% Preferred style for FDS and Smokeview plots 

paper_width     = 6.0;  % inches
paper_height    = 4.5;  % inches
plot_width      = 4.5;  % inches
plot_height     = 3.15; % inches
Key_Font_Size   = 12;
Title_Font_Size = 14;
Label_Font_Size = 14;
set(gca,'Units','inches')
set(gca,'FontName','Times')
set(gca,'Position',[1,0.75,plot_width,plot_height])