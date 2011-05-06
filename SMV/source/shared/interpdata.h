// $Date$ 
// $Revision$
// $Author$

/* --------------------------  xyzdata ------------------------------------ */

typedef struct {
  int *ijkbar;
  unsigned char *radiance, *opacity;
  float *xyzbar0, *xyzbar;
  float *dxyz;
} radiancedata;

typedef struct {
  float xyz[3];
} point;

typedef struct _kd_data {
  struct _kd_data *left, *right, *parent;
  int axis;
  point *median;
} kd_data;

typedef struct _scriptfiledata {
  struct _scriptfiledata *prev, *next;
  int id;
  int recording;
  char *file;
} scriptfiledata;


void setup_radiancemap(radiancedata *radianceinfo, int ijkbar[3], float xyzbar0[3], float xyzbar[3], float dxyz[3], 
                       unsigned char *radiance, unsigned char *opacity);
void build_radiancemap(radiancedata *radianceinfo);
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]);
void free_kdtree(kd_data *kdtree);
kd_data *setup_kdtree(point *points, int npoints, kd_data *parent);
