function[D] = vettori_sloped_function(file)
text=fileread(file);
for j =1:4
    name = ['TC_',num2str(j),'_1'];
    f1=findstr(name , text);
    D=textscan(text(f1(end)+30:f1(end)+40),'%[^ ]');
    format ('short', 'g');
    T(j,1)=str2num(D{1}{1});
end
format