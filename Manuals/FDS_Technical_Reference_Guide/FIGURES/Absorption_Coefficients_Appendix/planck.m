function I=planck(Temp,eta)
    if nargin<2
           eta=linspace(1,1.2e6,1000);
    end
    % Radiation constants
    C1=1.1910428220e-16; % W* m^2/sr
    C2=1.438775225e-2; % m*K
    I=((C1*eta.^3)./(exp((C2*eta)/Temp)-1)).';
end