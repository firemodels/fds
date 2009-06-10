% McDermott
% 6-02-2009
% dvcread.m
%
% Read a _devc.csv file for V&V processing.  The idea is to create the
% basic equivalent of csvread, except that we need also to read the header
% lines.  In principle, xlsread should handle this but we have experienced
% problems with this function.  Actually, the file given by 'filename' does
% not have to be a devc file.  But we take advantage of the fact that after
% the header lines have been read, the remaining data is numeric and 
% comma-delimited, so that we may use csvread at that point.  We also
% assume that the header names are contained in the last successfully read
% text line, since in devc files the units are typically given by the first
% line of text.

function [H X] = dvcread(filename)

fid = fopen(filename,'r+');

ierror = 0;
irow = -1;
textline = [];
while ierror==0
    irow = irow + 1;
    textline_last = textline;
    textline = fgetl(fid);
    [X ierror]=str2num(textline);
end
fclose(fid);

% parse header line
if irow>0
    C = textscan(textline_last,'%q','delimiter',',');
    H = strtrim(C{:}');
else
    H = [];
end  

% read numeric data
X = csvread(filename,irow,0);