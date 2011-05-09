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
  float dist2;
} kdpoint;

typedef struct _kd_data {
  struct _kd_data *left, *right, *parent;
  int axis;
  kdpoint *median;
} kd_data;

typedef struct _scriptfiledata {
  struct _scriptfiledata *prev, *next;
  int id;
  int recording;
  char *file;
} scriptfiledata;

kd_data *closest_node(kd_data *here, float *point, kd_data *best);
void sort_closest_nodes(kd_data **bests, int nbests, float *point);
float maxdistance2(kd_data **bests,int nbests, float *xyz,int *nmax);
kd_data *child_far(kd_data *here, float *point);
kd_data *child_near(kd_data *here, float *point);
int compare_pointz( const void *arg1, const void *arg2 );
float distance2(kd_data *node, float *xyz);
float distance_axis(kd_data *node, float *xyz);
int compare_pointy( const void *arg1, const void *arg2 );
int compare_pointy( const void *arg1, const void *arg2 );
int compare_pointx( const void *arg1, const void *arg2 );
int compare_points( const void *arg1, const void *arg2 );
int compare_bests( const void *arg1, const void *arg2 );
void build_radiancemap2(radiancedata *radianceinfo);
void test_kd(void);
void get_closest_points(kdpoint **pointers, int npoints, float *point);
void get_closest_nodes(kd_data *here, float *point, kd_data **bests, int *nbests, int nwanted);
void setup_radiancemap(radiancedata *radianceinfo, int ijkbar[3], float xyzbar0[3], float xyzbar[3], float dxyz[3], 
                       unsigned char *radiance, unsigned char *opacity);
void build_radiancemap(radiancedata *radianceinfo);
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]);
void free_kdtree(kd_data *kdtree);
kd_data *setup_kdtree(kdpoint *points, int npoints, kd_data *parent);
