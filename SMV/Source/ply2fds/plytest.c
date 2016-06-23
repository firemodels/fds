/*

Sample code showing how to read and write PLY polygon files.

Greg Turk, March 1994

*/

#include <stdio.h>
#include <math.h>
#include <ply.h>

/* user's vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;             /* the usual 3-space position of a vertex */
} Vertex;

typedef struct Face {
  unsigned char intensity; /* this user attaches intensity to faces */
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
} Face;

/* polygon description of an object (a cube) */

Vertex verts[] = {  /* vertices */
  { 0.0, 0.0, 0.0},
  { 1.0, 0.0, 0.0},
  { 1.0, 1.0, 0.0},
  { 0.0, 1.0, 0.0},
  { 0.0, 0.0, 1.0},
  { 1.0, 0.0, 1.0},
  { 1.0, 1.0, 1.0},
  { 0.0, 1.0, 1.0},
};

Face faces[] = {  /* faces */
  { '\001', 4, NULL },  /* intensity, vertex list count, vertex list (empty) */
  { '\004', 4, NULL },
  { '\010', 4, NULL },
  { '\020', 4, NULL },
  { '\144', 4, NULL },
  { '\377', 4, NULL },
};

/* list of vertices for each face */
/* (notice that indices begin at zero) */

typedef int Vertex_Indices[4];
Vertex_Indices vert_ptrs[] = {
  { 0, 1, 2, 3 },
  { 7, 6, 5, 4 },
  { 0, 4, 5, 1 },
  { 1, 5, 6, 2 },
  { 2, 6, 7, 3 },
  { 3, 7, 4, 0 },
};

/* information needed to describe the user's data to the PLY routines */

char *elem_names[] = { /* list of the kinds of elements in the user's object */
  "vertex", "face"
};

PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,z), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a vertex */
  {"intensity", PLY_UCHAR, PLY_UCHAR, offsetof(Face,intensity), 0, 0, 0, 0},
  {"vertex_indices", PLY_INT, PLY_INT, offsetof(Face,verts),
   1, PLY_UCHAR, PLY_UCHAR, offsetof(Face,nverts)},
};



/******************************************************************************
The main routine just creates and then reads a PLY file.
******************************************************************************/

main()
{

#if 1
  /* write a PLY file */
  write_test();
#endif

#if 1
  /* read a PLY file */
  read_test();
#endif

}


/******************************************************************************
Write out a PLY file.
******************************************************************************/

write_test()
{
  int i,j;
  PlyFile *ply;
  int nelems;
  char **elist;
  int file_type;
  float version;
  int nverts = sizeof (verts) / sizeof (Vertex);
  int nfaces = sizeof (faces) / sizeof (Face);

  /* create the vertex index lists for the faces */
  for (i = 0; i < nfaces; i++)
    faces[i].verts = vert_ptrs[i];

  /* open either a binary or ascii PLY file for writing */
  /* (the file will be called "test.ply" because the routines */
  /*  enforce the .ply filename extension) */

#if 1
  ply = ply_open_for_writing("test", 2, elem_names, PLY_ASCII, &version);
#else
  ply = ply_open_for_writing("test", 2, elem_names, PLY_BINARY_BE, &version);
#endif

  /* describe what properties go into the vertex and face elements */

  ply_element_count (ply, "vertex", nverts);
  ply_describe_property (ply, "vertex", &vert_props[0]);
  ply_describe_property (ply, "vertex", &vert_props[1]);
  ply_describe_property (ply, "vertex", &vert_props[2]);

  ply_element_count (ply, "face", nfaces);
  ply_describe_property (ply, "face", &face_props[0]);
  ply_describe_property (ply, "face", &face_props[1]);

  /* write a comment and an object information field */
  ply_put_comment (ply, "author: Greg Turk");
  ply_put_obj_info (ply, "random information");

  /* we have described exactly what we will put in the file, so */
  /* we are now done with the header info */
  ply_header_complete (ply);

  /* set up and write the vertex elements */
  ply_put_element_setup (ply, "vertex");
  for (i = 0; i < nverts; i++)
    ply_put_element (ply, (void *) &verts[i]);

  /* set up and write the face elements */
  ply_put_element_setup (ply, "face");
  for (i = 0; i < nfaces; i++)
    ply_put_element (ply, (void *) &faces[i]);

  /* close the PLY file */
  ply_close (ply);
}


/******************************************************************************
Read in a PLY file.
******************************************************************************/

read_test()
{
  int i,j,k;
  PlyFile *ply;
  int nelems;
  char **elist;
  int file_type;
  float version;
  int nprops;
  int num_elems;
  PlyProperty **plist;
  Vertex **vlist;
  Face **flist;
  char *elem_name;
  int num_comments;
  char **comments;
  int num_obj_info;
  char **obj_info;

  /* open a PLY file for reading */
  ply = ply_open_for_reading("test", &nelems, &elist, &file_type, &version);

  /* print what we found out about the file */
  printf ("version %f\n", version);
  printf ("type %d\n", file_type);

  /* go through each kind of element that we learned is in the file */
  /* and read them */

  for (i = 0; i < nelems; i++) {

    /* get the description of the first element */
    elem_name = elist[i];
    plist = ply_get_element_description (ply, elem_name, &num_elems, &nprops);

    /* print the name of the element, for debugging */
    printf ("element %s %d\n", elem_name, num_elems);

    /* if we're on vertex elements, read them in */
    if (equal_strings ("vertex", elem_name)) {

      /* create a vertex list to hold all the vertices */
      vlist = (Vertex **) malloc (sizeof (Vertex *) * num_elems);

      /* set up for getting vertex elements */

      ply_get_property (ply, elem_name, &vert_props[0]);
      ply_get_property (ply, elem_name, &vert_props[1]);
      ply_get_property (ply, elem_name, &vert_props[2]);

      /* grab all the vertex elements */
      for (j = 0; j < num_elems; j++) {

        /* grab and element from the file */
        vlist[j] = (Vertex *) malloc (sizeof (Vertex));
        ply_get_element (ply, (void *) vlist[j]);

        /* print out vertex x,y,z for debugging */
        printf ("vertex: %g %g %g\n", vlist[j]->x, vlist[j]->y, vlist[j]->z);
      }
    }

    /* if we're on face elements, read them in */
    if (equal_strings ("face", elem_name)) {

      /* create a list to hold all the face elements */
      flist = (Face **) malloc (sizeof (Face *) * num_elems);

      /* set up for getting face elements */

      ply_get_property (ply, elem_name, &face_props[0]);
      ply_get_property (ply, elem_name, &face_props[1]);

      /* grab all the face elements */
      for (j = 0; j < num_elems; j++) {

        /* grab and element from the file */
        flist[j] = (Face *) malloc (sizeof (Face));
        ply_get_element (ply, (void *) flist[j]);

        /* print out face info, for debugging */
        printf ("face: %d, list = ", flist[j]->intensity);
        for (k = 0; k < flist[j]->nverts; k++)
          printf ("%d ", flist[j]->verts[k]);
        printf ("\n");
      }
    }
    
    /* print out the properties we got, for debugging */
    for (j = 0; j < nprops; j++)
      printf ("property %s\n", plist[j]->name);
  }

  /* grab and print out the comments in the file */
  comments = ply_get_comments (ply, &num_comments);
  for (i = 0; i < num_comments; i++)
    printf ("comment = '%s'\n", comments[i]);

  /* grab and print out the object information */
  obj_info = ply_get_obj_info (ply, &num_obj_info);
  for (i = 0; i < num_obj_info; i++)
    printf ("obj_info = '%s'\n", obj_info[i]);

  /* close the PLY file */
  ply_close (ply);
}

