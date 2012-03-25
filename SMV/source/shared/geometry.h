// $Date$ 
// $Revision$
// $Author$

// svn revision character string
#define GEOM_INTERIOR 0
#define GEOM_EXTERIOR 1
#define GEOM_COMPLEX 2

/* --------------------------  vertdata ------------------------------------ */

typedef struct _vertdata {
  float xyz[3];
  int nedges,nverts,ntris;
  int type;
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
  vertdata *verts[3];
  int examined;
  float area,center[3],norm[3];
} tridata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _geomdata {
  vertdata *verts;
  tridata *tris;
  edgedata *edges;
  int nverts, ntris, nedges;
} geomdata;

