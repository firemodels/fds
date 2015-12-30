% McDermott
% 12-30-2015
% parseplus.m
%
% [S] = parseplus(character_string)
%
% Uses textscan to parse a character string that is delimited with a "+".

function [S] = parseplus(character_string)

cell_array = textscan(char(character_string),'%s','delimiter','+');
S = cell_array{:}';