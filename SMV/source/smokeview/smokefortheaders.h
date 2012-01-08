// $Date$ 
// $Revision$
// $Author$

#ifndef SMOKEFORTHEADERS_H_DEFINED
#define SMOKEFORTHEADERS_H_DEFINED

#include "isodefs.h"
#include "flowfiles.h"
#include "egz_stdio.h"

#ifdef WIN32
#define STDCALL extern void _stdcall
#else
#define STDCALL extern void
#endif
#ifdef pp_noappend
#define _F(name) name
#else
#define _F(name) name ## _
#endif

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
#define FORTgetdata2 _F(getdata2)
#define FORTgetsizes _F(getsizes)
#define FORTgetsizesa _F(getsizesa)
#define FORTgetsizes2 _F(getsizes2)
#define FORTgetslicesizes _F(getslicesizes)
#define FORTgetslicedata _F(getslicedata)
#define FORTgetplot3dq _F(getplot3dq)
#define FORTgetsliceparms _F(getsliceparms)
#define FORTcolor2rgb _F(color2rgb)
#define FORTget_file_unit _F(get_file_unit)
#define FORTclosefortranfile _F(closefortranfile) 
#define FORTgetboundaryheader1 _F(getboundaryheader1)
#define FORTgetboundaryheader2 _F(getboundaryheader2)


STDCALL FORTgeomout(float *verts, int *nverts, int *faces, int *nfaces);
STDCALL FORTgetembeddatasize(char *filename, int *endian, int *ntimes, int *nvars, int *error, FILE_SIZE lenfile);
STDCALL FORTgetembeddata(char *filename, int *endian, int *ntimes, int *nvals, float *times, int *nstatics, int *ndynamics,
                         float *vals, int *error, FILE_SIZE lenfile);
STDCALL FORTgetboundaryheader1(char *boundaryfilename, int *boundaryunitnumber, 
                               int *endian, int *npatch, int *error, FILE_SIZE lenfile);
STDCALL FORTgetboundaryheader2(int *boundaryunitnumber, int *version, int *npatches,
                               int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2, int *patchdir);
STDCALL FORTopenboundary(char *boundaryfilename, int *boundaryunitnumber, 
                         int *endian, int *version, int *error, FILE_SIZE len);
STDCALL FORTclosefortranfile(int *lunit);
STDCALL FORTget_file_unit(int *funit, int *first_unit);
STDCALL FORTcolor2rgb(int *rgb, char *color, FILE_SIZE colorsize);

STDCALL FORTfcreate_part5sizefile(char *part5file, char *part5sizefile, int *angle_flag, int *error,
                                  FILE_SIZE lenpart5file, FILE_SIZE lenpart5sizefile);

STDCALL FORTgetsliceparms(char *file,int *endian,
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,int *ni, int *nj, int *nk, int *slice3d, int *error,FILE_SIZE lenfile);
STDCALL FORTgetzonesize(char *zonefilename, int *nzonet, int *nrooms, int *nfires, int *endian, int *error, FILE_SIZE len);
STDCALL FORTgetzonedata(char *zonefilename, int *nzonet, int *nrooms, int *nfires, 
                        float *zonet, float *zoneqfire, float *zonepr, float *zoneylay,float *zonetl,float *zonetu, int *endian,
                        int *error, FILE_SIZE len);
STDCALL FORTgetxyzdata(int *iblank,int *nx,int *ny,int *nz,int *error);
STDCALL FORTgetpatchsizes1(int *file_unit,char *patchfilename,char *patchlonglabel,char *patchshortlabel,char *patchunit,
                           int *endian,int *npatch, int *headersize, int *error,
                           FILE_SIZE len1, FILE_SIZE len2, FILE_SIZE len3, FILE_SIZE len4);
STDCALL FORTgetpatchsizes2(int *file_unit,int *version, int *npatch,int *npatchsize, 
                           int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2, int *patchdir,
                           int *headersize, int *framesize);
STDCALL FORTgetpatchdata(int *lunit, int *npatch,int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2,
                         float *patchtimes,float *pqq, int *npqq, int *error);
STDCALL FORTgetdata1(int *file_unit, int *ipart, int *error);
STDCALL FORTgetdata2(int *file_unit,
                     short *xparts, short *yparts, short *zparts,
                     float *t,
                     int *sprinkflag, unsigned char *isprink, float *tspr, int *bframe,int *sframe,int *sprframe,
                     float *ptimes,int *nspr,int *nmax,int *mxframes,int *nframes,
                     int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p, int *partframestep, int *partpointstep, 
                     float *xbox0, float *xbox, float *ybox0, float *ybox, float *zbox0, float *zbox, 
                     float *offset_x, float *offset_y, float *offset_z,
                     int *error,
                     FILE_SIZE lenisprink);

STDCALL FORTgetsizesa(char *partfilename, int *npartpoint,int *npartframes,FILE_SIZE lenfile);
STDCALL FORTgetsizes(int *file_unit,char *partfilename, int *ibar, int *jbar, int *kbar, 
                     int *nb, int *nv, int *nspr,int *mxframepoints, int *endian, int *showstaticsmoke, int *error, FILE_SIZE filelen);
STDCALL FORTgetsizes2(int *file_unit,int *settmin_p, float *tmin_p, int *settmax_p, float *tmax_p,
                      int *nspr, int *frameloadstep, int *partpointstep, int *npartpoints, int *npartframes, int *error);
STDCALL FORTgetslicesizes(char *slicefilename, int *nslicei, int *nslicej, int *nslicek, 
                          int *nsteps,int *sliceframestep, int *endian, int *error, 
                          int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p,
                          int *headersize, int *framesize, int *statfile,
                          FILE_SIZE slicefilelen);
STDCALL FORTgetslicedata(int *file_unit,char *slicefilename, char *longlabels, char *shortlabels, 
                         char *units, int *is1,int *is2,int *js1,int *js2, int *ks1,int *ks2,
                         int *idir, float *qslicemin,float *qslicemax,
                         float *qslicedata,float *slicetimes, int *nsteps, int *sliceframestep, int *endian,
                         int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p,
                         FILE_SIZE slicefilelen, FILE_SIZE longlableslen, FILE_SIZE shortlabelslen, FILE_SIZE unitslen);
STDCALL FORTgetplot3dq(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, int *endian, int *isotest, FILE_SIZE filelen);
#endif
