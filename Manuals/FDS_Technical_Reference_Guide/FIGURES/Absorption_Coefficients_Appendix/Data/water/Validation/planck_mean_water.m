waterk = readdata('water_k.txt',1);
lambda = 1:200;
T = 300:20:1000;
for t = 1:length(T)
   for l = 1:length(lambda)
      n_i = interp1(waterk(:,1),waterk(:,2),lambda(l)/1E6);
      kw(l) = 4*pi*n_i/(lambda(l)*1E-6);
      I(l) = planck(lambda(l),T(t));
   end
   k(t)=trapz(lambda,I.*kw)/(5.67E-8*T(t)^4/pi);
   icheck(t)=trapz(lambda,I)/(5.67E-8*T(t)^4/pi);
end
Estr = 'Mean relative error of Planck integration: ';
disp([Estr num2str(mean(abs(1-icheck)))])
