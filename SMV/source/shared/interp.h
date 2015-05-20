#ifndef INTERP_H_DEFINED
#define INTERP_H_DEFINED
#define IJK(i,j,k) ((i)+(j)*nx+(k)*nxy)
void get_z_interp_factors(float *zplt, int nz, float z, int *k1, int *k2, float *f1, float *f2);
int interp3dsliceindex(unsigned char *data, float *zplt, int nz, int n0, float z);
#endif
