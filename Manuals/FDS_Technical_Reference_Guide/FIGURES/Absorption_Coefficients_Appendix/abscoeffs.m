
function [a v]=abscoeffs(data)
if(isstruct(data))
    if(isfield(data,'type'))
        if(strcmp(data.type,'jcamp_struct'))
            D=data;
        else
            error('Wrong type struct');
        end
    else
        error('Wrong type struct');
    end
else
    if(ischar(data))
        D=readjcamp(data);
    else
        error('First argument needs to be jcamp_struct or filename');
    end
end


% 1/cm=1/(m*10^-2)=1/m*10^2
v=D.x*10^2;  % wave number 1/cm->1/m


sigma=5.67040040*10^-8; 
S=D.path_length*1e-2; % path length (cm=>m)
foo=regexp(D.state,'\w*','match');
state=foo{1};
switch(upper(state))
    case 'SOLUTION'
        c=str2double(foo(2))*10; % wt/vol% (g/100ml => 10 kg/m^3)
        aunits='$m^3 \cdot kg^{-1} \cdot m^{-1}$';
    case 'GAS'
        if(length(foo)>1)
            c=str2double(foo(2))*133.322; %mmHG -> Pa
            aunits='$(Pa\cdot m)^{-1}$';
        else
            c=1;
            aunits='$ m^{-1}$';
        end
    case 'LIQUID'
        c=1;
        aunits='$ m^{-1}$';
    otherwise
        error(strcat('Unknown state: ',foo(1)));
end



switch(upper(D.yunits))
    case 'TRANSMITTANCE'

        T=D.y; % Transmittance
        %T(T<=0)=eps;
        T(T<=0)=realmin;
        A=-log(T);
    case 'ABSORBANCE'
        
        A=D.y; %Absorbance
end
%Spectral absorption coeff.
a=A./(S*c);
a(a<0)=0;
end