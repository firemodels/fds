% McDermott
% 6-08-2009
% parse.m
%
% [S] = parse(character_string)
%
% Uses textscan to parse a character string that is delimited with a "|".

function [S] = parse(character_string)

cell_array = textscan(character_string,'%s','delimiter','|');
S = cell_array{:}';