% McDermott
% 7-21-11
% wavelet_transform.m
%
% Matlab version of WAVELET_ERROR function in dump.f90

close all
clear all

S = [0 1 0 1];

plot([1:length(S)],S,'bo--')
xlabel('cell index')
ylabel('normalized signal')

% normalize signal
SMAX=max(S);
SMIN=min(S);
DS=SMAX-SMIN;
if DS<1.E-6
    WAVELET_ERROR = 0.;
else
    SS=(S-SMIN)/DS;
    
    % discrete Haar wavelet transform
    M=2;
    N=M;
    for I=1:M
        for J=1:N
            K=2*J-1;
            if I==1
                A(I,J) = 0.5*(SS(K)+SS(K+1));
                C(I,J) = 0.5*(SS(K)-SS(K+1));
            else
                A(I,J) = 0.5*(A(I-1,K)+A(I-1,K+1));
                C(I,J) = 0.5*(A(I-1,K)-A(I-1,K+1));
            end
        end
        N=N/2;
    end
    
    C1 = sum(C(1,:));
    C2 = sum(C(2,:));
    
    WAVELET_ERROR = abs(C1-C2)
    A
    C
end

