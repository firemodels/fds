% McDermott
% 6-08-2009
% parsepipe.m
%
% [S] = parsepipe(character_string)
%
% Uses textscan to parse a character string that is delimited with a "|".

function [S] = parsepipe(character_string)

cell_array = textscan(character_string,'%s','delimiter','|');
S = cell_array{:}';