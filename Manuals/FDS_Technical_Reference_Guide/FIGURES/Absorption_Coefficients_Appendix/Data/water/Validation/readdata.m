function M = readdata(FNAME,NI,NROWS, NCOLS)
% READDATA      Read data from ascii file
%
% M = READDATA(FNAME,Ni[,NROWS,NCOLS])
%
% reads data to matrix M from file FNAME ignoring first NI rows
% Ni is the number of header rows.
% If the number of data columns changes, stops reading.
%
if (nargin<3) 
   NROWS = 1;
   NCOLS = 1;
end
fid = fopen(FNAME,'r');
if (NI>0)
for i = 1:NI
   whorow = fgetl(fid);
end
end
if (NROWS<=1)
   M = [];
else
   M(NROWS,NCOLS) = 0;
end
row = fgetl(fid);
i = 0;
while row ~= -1
    i = i + 1;
    if (eq(mod(i,10000),0))
       disp(i)
    end	
%   if row(1) == 'S'
%      row = row(3:length(row)-3);
%   end
   MM = str2num(str2mat(row));
   if (i>1)
      if (length(MM)~=length(M(i-1,:))), break; end
   end
   M = [M; MM];
   row = fgetl(fid);
end
fclose(fid);
