// $Date$ 
// $Revision$
// $Author$

#ifndef KDTREE_H_DEFINED
#define KDTREE_H_DEFINED
#ifdef INKDTREE
#define KDEXTERN
#else
#define KDEXTERN extern
#endif

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

KDEXTERN kd_data *setup_kdtree(kdpoint *points, int npoints, kd_data *parent);
KDEXTERN void get_closest_nodes(kd_data *here, float *point, kd_data **bests, int *nbests, int nwanted);
KDEXTERN void sort_closest_nodes(kd_data **bests, int nbests, float *point);
KDEXTERN void get_closest_points(kdpoint **pointers, int npoints, float *point);
#endif
