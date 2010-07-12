function[T] = vettori_function(file)
%addpath('Y:\FDS-SMV\Validation');
text=fileread(file);

for j =1:4
    name = ['SPRINKLER ',num2str(j)];
    f1=findstr(name , text);
    C=textscan(text(f1(end)+30:f1(end)+40),'%[^ ]');
    format ('short', 'g');
    T(j,1)=str2num(C{1}{1});
end
format