function x0 = inverse_percentile(x,y)

n = length(x);

for i=2:n
   psum(i) = trapz(x(1:i),y(1:i));
end
psum(1) = 0.;

k = min(find(psum>0.90*psum(n)));
if (k>1)
   x0 = x(k-1) + (x(k)-x(k-1))*(0.90*psum(n)-psum(k-1))/(psum(k)-psum(k-1));
else
   x0 = x(1);
end

