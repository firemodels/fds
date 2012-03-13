// $Date: 2012-03-10 17:31:20 -0500 (Sat, 10 Mar 2012) $ 
// $Revision: 10233 $
// $Author: gforney $

// svn revision character string

/* --------------------------  vertdata ------------------------------------ */

typedef struct _vertdata {
  float xyz[3];
  int nvertring,ntris;
  struct _vertdata **vertring;
  struct _tridata **tris;
  float best_xyz[3], best_norm[3];
} vertdata;

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
  int nverts, ntris;
} geomdata;

