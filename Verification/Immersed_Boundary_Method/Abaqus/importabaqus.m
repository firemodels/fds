% McDermott
% 10-19-11
% importabaqus.m
%
% Read nodes and volumes from Abaqus input file
%
% Example: [N V] = importabaqus(filename)

function [N V] = importabaqus(filename)

if ~exist(filename)
    display('file does not exist')
    return
end

fid = fopen(filename);

% skip header info

while 1
    tline = fgetl(fid);
    %disp(tline)
    if strcmp(tline,'*Node'); break; end
end

% read node lines

i=0;
while 1
    tline = fgetl(fid);
    if strcmp(tline(1),'*'); break; end
    nodes = str2num(tline);
    i=i+1;
    N(i,1:3) = nodes(2:4);
end

% read element lines

i=0;
while 1
    tline = fgetl(fid);
    if strcmp(tline(1),'*'); break; end
    nodes = str2num(tline);
    i=i+1;
    V(i,1:4) = nodes(2:5);
end

fclose(fid);
