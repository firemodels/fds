function En = expintn(n,x)
% EXPINTN(n,x)
%       ============================================
%       Purpose: Compute exponential integral En(x)
%       Input :  x --- Argument of En(x,x ó 20)
%                n --- Order of En(x)
%       Output:  EN(n)--- En(x)
%       Routine called: E1XB for computing E1(x)
%       ============================================
en(1)=exp(-x)./x;
e1=e1xb(x);
en(1+1)=e1;
for  k=2:n;
   ek=(exp(-x)-x.*e1)./(k-1.0d0);
   en(k+1)=ek;
   e1=ek;
end;  
k=fix(n)+1;
En = en(n+1);
return;
end

function e1=e1xb(x)
%       ============================================
%       Purpose: Compute exponential integral E1(x)
%       Input :  x  --- Argument of E1(x)
%       Output:  E1 --- E1(x)
%       ============================================
if(x == 0.0)
   e1=1.0d+300;
elseif(x <= 1.0)
   e1=1.0;
   r=1.0;
   for  k=1:25
      r = -r.*k.*x./(k+1.0d0).^2;
      e1=e1+r;
      if(abs(r)<= (abs(e1).*1.0d-15)) 
         break; 
      end
   end
   ga=0.5772156649015328d0;
   e1=-ga-log(x)+x.*e1;
else
   m=20+fix(80.0./x);
   t0=0.0d0;
   for  k=m:-1:1;
      t0=k./(1.0d0+k./(x+t0));
   end
   t=1.0d0./(x+t0);
   e1=exp(-x).*t;
end
return
end


