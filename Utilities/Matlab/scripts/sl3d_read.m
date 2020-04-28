% McDermott
% 15 Aug 2019
% sl3d_read.m
%
% based on slread.m by Simo Hostikka
%
% SL3D_READ  Reads an FDS 3D slice file
% [Q,Time]=sl3d_read(fname,Tstart,Tend [,Tstep]);
%            Tstart  is start time
%            Tend    is end time
%            Tstep   is time step of dumps
%
%            Q       contains the data
%            Time    contains the time points
function [Q,Time]=sl3d_read(fname,Tstart,Tend,Tstep)

if nargin<4; Tstep = Tend-Tstart; end

Time = 0;

fid  = fopen(fname,'r','l'); % might need 'b' for "big endian" instead of 'l', check your file system

r4b(fid);
Str1 = char([fread(fid,30,'schar')]')
r4b(fid);
r4b(fid);
Str2 = char([fread(fid,30,'schar')]')
r4b(fid);
r4b(fid);
Str3 = char([fread(fid,30,'schar')]')
r4b(fid);

r4b(fid);
Indx = [fread(fid,6,'int32')]'
r4b(fid);

% allocate Q

Isize = Indx(2)-Indx(1)+1;
Jsize = Indx(4)-Indx(3)+1;
Ksize = Indx(6)-Indx(5)+1;

Nrun = max(1,round((Tend-Tstart)/Tstep));
Q(Isize,Jsize,Ksize,Nrun) = 0;
Q_tmp(Isize,Jsize*Ksize) = 0;

st = 1;

while Time < Tstart

    % disp(['reading Time ' num2str(st)])

    r4b(fid);
    Time(st) = fread(fid,1,'float32');
    r4b(fid);

    r4b(fid);
    Q_tmp = fread(fid,[Isize*Jsize*Ksize],'float32')';
    Q(:,:,:,st) = reshape(Q_tmp,[Isize,Jsize,Ksize]);
    r4b(fid);

end

while Time < Tend,

    st = st + 1;

    % disp(['reading Time ' num2str(st)])

    r4b(fid);
    Time(st) = fread(fid,1,'float32');
    r4b(fid);

    r4b(fid);
    Q_tmp = fread(fid,[Isize*Jsize*Ksize],'float32')';
    Q(:,:,:,st) = reshape(Q_tmp,[Isize,Jsize,Ksize]);
    r4b(fid);

end

fclose(fid);

function r4b(fid)
fread(fid,4,'int8');