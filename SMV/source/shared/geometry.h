// $Date: 2012-03-10 17:31:20 -0500 (Sat, 10 Mar 2012) $ 
// $Revision: 10233 $
// $Author: gforney $

// svn revision character string

/* --------------------------  vertdata ------------------------------------ */

typedef struct _vertdata {
  float xyz[3];
  int nvertring,ntrifan;
  struct _vertdata **vertring;
  struct _edgedata **edgering;
  struct _tridata **trifan;
  float best_xyz[3], best_norm[3];
} vertdata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _edgedata {
  vertdata *verts[2];
  struct _tridata **tris;
  int ntris;
  float center[3];
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

