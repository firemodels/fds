% Simo Hostikka
% 20 Oct 2008
% slread.m
%
% SLREAD    Reads a FDS slice file
% [T,Time]=slread(fname,Tstart,Tend [,Tstep]);
%           Tstart  is start time
%           Tend    is end time
%           Tstep   is time step of dumps
%
%       T       contains the data
%       Time    contains the time points
function [T,Time]=slread(fname,Tstart,Tend,Tstep)

if (nargin<4),Tstep = Tend-Tstart;end

Time = 0;

if (filesep=='/')
    fid  = fopen(fname,'r','b');
else
    fid  = fopen(fname,'r','l');
end
r4b(fid);
Str1 = char([fread(fid,30,'schar')]');
r4b(fid);
r4b(fid);
Str2 = char([fread(fid,30,'schar')]');
r4b(fid);
r4b(fid);
Str3 = char([fread(fid,30,'schar')]');
r4b(fid);

r4b(fid);
Indx = [fread(fid,6,'int32')]';
r4b(fid);

% allocate T

Isize = Indx(2)-Indx(1)+1;
Jsize = Indx(4)-Indx(3)+1;
Ksize = Indx(6)-Indx(5)+1;
if (Isize == 1)
   M = Jsize;
   N = Ksize;
elseif (Jsize == 1)
   M = Isize;
   N = Ksize;
else
   M = Isize;
   N = Jsize;
end

Nrun = max(1,round((Tend-Tstart)/Tstep));
T(N,M,Nrun) = 0;

st = 1;

while Time < Tstart,
r4b(fid);
Time(st) = fread(fid,1,'float32');
r4b(fid);

r4b(fid);
T(:,:,st) = fread(fid,[M,N],'float32')';
r4b(fid);

end

while Time < Tend,

st = st + 1;

r4b(fid);
Time(st) = fread(fid,1,'float32');
r4b(fid);

r4b(fid);
T(:,:,st) = fread(fid,[M,N],'float32')';
r4b(fid);

end

%T = [T/st]';
%if st > 1,
%  disp(['Number of averaging steps = ' num2str(st) '.'])
%end

fclose(fid);

function r4b(fid)
fread(fid,4,'int8'); 