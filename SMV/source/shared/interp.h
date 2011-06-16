// $Date$ 
// $Revision$
// $Author$


#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]);
