% McDermott
% 6-11-2009
% stripcell.m
%
% Stip cell array of empty cells

function [C] = stripcell(CELL_ARRAY)

i = 1;
for j=1:length(CELL_ARRAY)
    if iscellstr(CELL_ARRAY(j))
        C(i) = CELL_ARRAY(j);
        i = i+1;
    end
end