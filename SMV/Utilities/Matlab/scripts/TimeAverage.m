%T. Myers
%6-07-2010 
%(Modified by A. Matala 07-09-2013)
%TimeAverage.m
%
%Input: in = file_name.csv; Intended to be used with a csv format file. The
%first x rows can be title or heading, which will be reinserted at the top
%of the output file. Cannot handle text down an entire column or text interspaced with the data. All text must be confined to title/heading rows
%
%Averages every f data points and creates a new data point at the end of
%the time period examined. Assumes time is in seconds.
%
%Output: out = file_name.csv; will out put a csv format file. The exact
%title and heading rows will be inserted into the initial rows, and the
%adjusted data will be appended to the end.


clear, clc;
%Determine the file from which to draw the data
in = input('Please enter the name of the file to be time averaged. \n Be sure the file is contained in MATLAB’s current folder. \n :','s');

%Pulls data from the given file, beginning in given row
x = input('Which row does the data begin in? \n :');
M=csvread(in,(x-1),0);

%Requests the factor by which data will be time requested. Error will be
%produced if the number of data steps is not evenly divisible by the
%factor.
f = input('Please specify the factor by which you would like to reduce the data. \n :');

%Sizes the data file, and creates a new matrix to populate with the time
%averaged data
[v,u]=size(M);

%Interpolate to make dt = 1 s
%First, make sure that the values of time column are distinct
d_indx = find(diff(M(:,1))>0);

M2(:,1) = [ceil(M(d_indx(1),1)):floor(M(d_indx(length(d_indx)),1))];
for i = 2:u
   M2(:,i) = interp1(M(d_indx,1),M(d_indx,i),M2(:,1));
end

N=zeros(floor(length(M2(:,1))/f)-1,u);
%Time averages the data, first creating the time row, then creating a zero
%point, then running the averages.
N(1,:) = M2(1,:);
t_min = N(1,1);
for m =2:u
for n =1:length(N(:,1))
%     N(1,1)=0;
%     N(n+1,1)=((f)*(n));
%     N(1,m)=M2(1,m);
%     N(n+1,m)=mean(M2((1+(f)*(n-1)):((f)*n),m));
    N(n,1) = t_min + f*n;
    
    if isequal(rem(f,2),0) %even numbers
        N(n,m) = mean(M2((f*n-f/2+1):(f*n+f/2),m));
    
    else % odd numbers
        N(n,m) = mean(M2((f*n-floor(f/2)):(f*n+floor(f/2)),m));
    end
end
end

%round the results to requested presision
p = input('Please specify the presision of the results.\n If the columns have different presision, give them in brackets []. \n :');
if isequal(length(p),1)
    N = round(N./p).*p;
else
    for i = 1:length(p)
       N(:,i) =  round(N(:,i)./p(i)).*p(i);
    end
end

%Re-reads the data file, this time to acquire the text data. 
[A,B,J] = xlsread(in);

%Requests a name for an output file. Regardless of file type specified a
%file in csv format will be generated.
out = input('Please enter the desired output file name \n The file type must be included in the name \n :','s');

%Writes the numerical data, leaving a space for text. Then attaches
%the text in the space left
dlmwrite(out,N,',',(x-1),0);
xlswrite(out,B);

%Confirms the data has been stored in given location
disp( ['Time averaged data is stored in ' out] )


