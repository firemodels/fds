#ifndef SMOKEFORTHEADERS_H_DEFINED
#define SMOKEFORTHEADERS_H_DEFINED

#ifdef WIN32
#define STDCALLF extern void _stdcall
#else
#define STDCALLF extern void
#endif

#define FORTfpoly2tri _F(fpoly2tri)
#define FORTget_in_triangle _F(get_in_triangle)
#define FORTget_is_angle_ge_180 _F(get_is_angle_ge_180)
#define FORTtest_in_tetra _F(test_in_tetra)
#define FORTgetverts _F(getverts2)
#define FORTgettetravol _F(get_tetrabox_volume_fb)
#define FORTgeomout _F(geomout)
#define FORTgetembeddatasize _F(getembeddatasize)
#define FORTgetembeddata _F(getembeddata)
#define FORTopenboundary _F(openboundary)
#define FORTfcreate_part5sizefile _F(fcreate_part5sizefile)
#define FORTgetzonesize _F(getzonesize)
#define FORTgetzonedata _F(getzonedata)
#define FORTgetxyzdata _F(getxyzdata)
#define FORTgetpatchsizes1 _F(getpatchsizes1)
#define FORTgetpatchsizes2 _F(getpatchsizes2)
#define FORTgetpatchdata _F(getpatchdata)
#define FORTgetdata1 _F(getdata1)
#define FORTgetsizes _F(getsizes)
#define FORTgetsizesa _F(getsizesa)
#define FORTgetslicesizes _F(getslicesizes)
#define FORTwriteslicedata _F(writeslicedata)
#define FORTgetslicedata _F(getslicedata)
#define FORTgetplot3dq _F(getplot3dq)
#define FORTgetsliceparms _F(getsliceparms)
#define FORTcolor2rgb _F(color2rgb)
#define FORTget_file_unit _F(get_file_unit)
#define FORTclosefortranfile _F(closefortranfile)
#define FORTgetboundaryheader1 _F(getboundaryheader1)
#define FORTgetboundaryheader2 _F(getboundaryheader2)

STDCALLF FORTfpoly2tri(float *verts, int *nverts, int *poly, int *npoly, int *tris, int *ntris);
STDCALLF FORTget_in_triangle(float *vert, float *v1, float *v2, float *v3, int *flag);
STDCALLF FORTget_is_angle_ge_180(float *v1, float *v2, float *v3, int *flag);
STDCALLF FORTtest_in_tetra(float *xyz, int *in_tetra, int *tetra_state);
STDCALLF FORTgettetravol(float *box_bounds,float *v0,float *v1,float *v2,float *v3,float *tetra_vol,float *areas,float *centroid);
STDCALLF FORTgetverts(float *box_bounds, float *v0, float *v1, float *v2, float *v3, float *out_verts,
                      int *nverts, int *faces, int *face_id, int *which_poly, int *nfaces, int *npolys, int *box_state);
STDCALLF FORTgeomout(float *verts, int *nverts, int *faces, int *nfaces);
STDCALLF FORTgetembeddatasize(char *filename, int *ntimes, int *nvars, int *error, FILE_SIZE lenfile);
STDCALLF FORTgetembeddata(char *filename, int *ntimes, int *nvals, float *times, int *nstatics, int *ndynamics,
                         float *vals, int *redirect, int *error, FILE_SIZE lenfile);
STDCALLF FORTgetboundaryheader1(char *boundaryfilename, int *boundaryunitnumber,
                               int *npatch, int *error, FILE_SIZE lenfile);
STDCALLF FORTgetboundaryheader2(int *boundaryunitnumber, int *version, int *npatches,
                               int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2, int *patchdir);
STDCALLF FORTopenboundary(char *boundaryfilename, int *boundaryunitnumber,
                         int *version, int *error, FILE_SIZE len);
STDCALLF FORTclosefortranfile(int *lunit);
STDCALLF FORTget_file_unit(int *funit, int *first_unit);
STDCALLF FORTcolor2rgb(int *rgb, char *color, FILE_SIZE colorsize);

STDCALLF FORTfcreate_part5sizefile(char *part5file, char *part5sizefile, int *angle_flag, int *redirect_flag, int *error,
                                  FILE_SIZE lenpart5file, FILE_SIZE lenpart5sizefile);

STDCALLF FORTgetsliceparms(char *file,
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,int *ni, int *nj, int *nk, int *slice3d, int *error,FILE_SIZE lenfile);
STDCALLF FORTgetzonesize(char *zonefilename, int *nzone_times, int *nrooms, int *nfires, int *error, FILE_SIZE len);
STDCALLF FORTgetzonedata(char *zonefilename, int *nzone_times, int *nrooms, int *nfires,
                        float *zone_times, float *zoneqfire, float *zonepr, float *zoneylay,float *zonetl,float *zonetu,
                        int *error, FILE_SIZE len);
STDCALLF FORTgetxyzdata(int *iblank,int *nx,int *ny,int *nz,int *error);
STDCALLF FORTgetpatchsizes1(int *file_unit,char *patchfilename,char *patchlonglabel,char *patchshortlabel,char *patchunit,
                           int *npatch, int *headersize, int *error,
                           FILE_SIZE len1, FILE_SIZE len2, FILE_SIZE len3, FILE_SIZE len4);
STDCALLF FORTgetpatchsizes2(int *file_unit,int *version, int *npatch,int *npatchsize,
                           int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2, int *patchdir,
                           int *headersize, int *framesize);
STDCALLF FORTgetpatchdata(int *lunit, int *npatch,int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2,
                         float *patch_times,float *pqq, int *npqq, int *error);
STDCALLF FORTgetdata1(int *file_unit, int *ipart, int *error);

STDCALLF FORTgetsizesa(char *partfilename, int *npartpoint,int *npartframes,FILE_SIZE lenfile);
STDCALLF FORTgetsizes(int *file_unit,char *partfilename,
                     int *nb, int *nv, int *nspr,int *mxframepoints, int *showstaticsmoke, int *error, FILE_SIZE filelen);
STDCALLF FORTgetslicesizes(char *slicefilename, int *nslicei, int *nslicej, int *nslicek,
                          int *nsteps,int *sliceframestep, int *error,
                          int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p,
                          int *headersize, int *framesize,FILE_SIZE slicefilelen);
STDCALLF FORTwriteslicedata(int *file_unit,char *slicefilename,
                            int *is1,int *is2,int *js1,int *js2,int *ks1,int *ks2,
                            float *qdata,float *times,int *ntimes, int *redirect,FILE_SIZE slicefilelen);
STDCALLF FORTgetslicedata(int *file_unit,char *slicefilename, char *shortlabels,
                          int *is1,int *is2,int *js1,int *js2, int *ks1,int *ks2,
                          int *idir, float *qslicemin,float *qslicemax,
                          float *qslicedata,float *times, int *nsteps, int *sliceframestep,
                          int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p, int *redirect,
                          FILE_SIZE slicefilelen, FILE_SIZE shortlabelslen);
STDCALLF FORTgetplot3dq(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, int *isotest, FILE_SIZE filelen);
#endif
