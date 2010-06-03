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
N=zeros(((v./f)+1),u);

%Time averages the data, first creating the time row, then creating a zero
%point, then running the averages.
for m =2:u
for n =1:(v./f)
    N(1,1)=0;
    N(n+1,1)=((f)*(n));
    N(1,m)=M(1,m);
    N(n+1,m)=mean(M((1+(f)*(n-1)):((f)*n),m));
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


