#define GEOM_INTERIOR 0
#define GEOM_EXTERIOR 1
#define GEOM_COMPLEX 2

/* --------------------------  vertdata ------------------------------------ */

typedef struct _vertdata {
  float xyz[3];
  int nedges,nverts,ntris;
  int type,active,seq_id;
  struct _vertdata *prev,*next;
  struct _vertdata **verts;
  struct _edgedata **edges;
  struct _tridata **tris;
} vertdata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _edgedata {
  vertdata *verts[2];
  int ntris;
} edgedata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _tridata {
  struct _tridata *prev,*next;
  vertdata *verts[3];
  int examined,active,seq_id;
  float area,center[3],norm[3];
} tridata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _geomdata {
  vertdata *verts,*first_vert,*last_vert;
  tridata *tris,*first_tri,*last_tri;
  edgedata *edges;
  int nverts, ntris, nedges;
} geomdata;

int in_sphere(float *pt, float *center, float radius);
int in_cylinder(float *pt, float *base, float h, float radius);
int in_cone(float *pt, float *base, float h, float radius);


