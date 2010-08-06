function[T] = vettori_sloped_function(file)
text=fileread(file);
for j =1:4
    name = ['SPK_',num2str(j)];
    f1=findstr(name , text);
    C=textscan(text(f1(end)+25:f1(end)+40),'%[^ ]');
    format ('short', 'g');
    T(j,1)=str2num(C{1}{1});
end
format