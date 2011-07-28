cd U:\FISU\Tutkimus\3_Palokuormat\AbsorptionCoefficients\Data\water\Validation
planck_mean_water
cd U:\FISU\Tutkimus\3_Palokuormat\AbsorptionCoefficients\Scripts
foo = readjcamp('..\Data\Water\7732-18-5-IR.jdx');
[plank]=planckmean(foo,T);
cd U:\FISU\Tutkimus\3_Palokuormat\AbsorptionCoefficients\Data\water\Validation
h=plot(T,plank,T,k);
figformat(h,gca,gcf,'T (K)','kappa (m^{-1})','H_2O',0.8,12,2)
legend('NIST data','Hale & Querry')
figsaver('H2O_NIST_vs_HaleQuerry')
