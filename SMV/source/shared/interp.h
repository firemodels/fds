// $Date$ 
// $Revision$
// $Author$

#define IJK(i,j,k) ((i)+(j)*nx+(k)*nxy)
#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
float interp3d(float xyz[3], float *vals, mesh *meshi, int *inobst, char *blank);
void get_z_interp_factors(float *zplt, int nz, float z, int *k1, int *k2, float *f1, float *f2);
int interp3dsliceindex(unsigned char *data, float *zplt, int nz, int n0, float z);
