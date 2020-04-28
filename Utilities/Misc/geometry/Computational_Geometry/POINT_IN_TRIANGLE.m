function [PTINT_FLG]=POINT_IN_TRIANGLE(P,V1,V2,V3)

EPS=1.E-15;
PTINT_FLG=true;

% compute face normal
E1 = V2-V1;
E2 = V3-V1;
[N]=cross(E1,E2);

for I=1:3
    switch(I)
        case(1)
            E = V2-V1;
            R = P-V1;
        case(2)
            E = V3-V2;
            R = P-V2;
        case(3)
            E = V1-V3;
            R = P-V3;
    end
    Q=cross(E,R);
    if ( dot(Q,N) < -EPS )
        PTINT_FLG=false;
        return
    end
end


return