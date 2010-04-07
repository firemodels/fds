// $Date$ 
// $Revision$
// $Author$

#ifndef DEF_smokefortheaders
#define DEF_smokefortheaders

#include "isodefs.h"
#include "flowfiles.h"
#include "egz_stdio.h"

#ifdef WIN32
#define STDCALL extern void _stdcall
#else
#define STDCALL extern void
#endif
#ifndef pp_noappend
#define FORTfcreate_part5sizefile fcreate_part5sizefile_
#define FORTclosepart closepart_
#define FORTclosezone closezone_
#define FORTclosepatch closepatch_
#define FORTgetzonesize getzonesize_
#define FORTgetzonedata getzonedata_
#define FORTgetxyzsize getxyzsize_
#define FORTgetxyzdata getxyzdata_
#define FORTgetpatchsizes1 getpatchsizes1_
#define FORTgetpatchsizes2 getpatchsizes2_
#define FORTgetpatchdata getpatchdata_
#define FORTgetdata1 getdata1_
#define FORTgetdata2 getdata2_
#define FORTgetdata2b getdata2b_
#define FORTgetdata3 getdata3_
#define FORTgetdata2a getdata2a_
#define FORTgetsizes getsizes_
#define FORTgetsizes2 getsizes2_
#define FORTgetsizesa getsizesa_
#define FORTgetslicesizes getslicesizes_
#define FORTgetslicedata getslicedata_
#define FORTgetplot3dq getplot3dq_
#define FORTgetplot3dqa getplot3dqa_
#define FORTgetsliceparms getsliceparms_
#define FORTcolor2rgb color2rgb_
#else
#define FORTfcreate_part5sizefile fcreate_part5sizefile
#define FORTclosepart  closepart
#define FORTclosezone closezone
#define FORTclosepatch closepatch
#define FORTgetzonesize getzonesize
#define FORTgetzonedata getzonedata
#define FORTgetxyzsize getxyzsize
#define FORTgetxyzdata getxyzdata
#define FORTgetpatchsizes1 getpatchsizes1
#define FORTgetpatchsizes2 getpatchsizes2
#define FORTgetpatchdata getpatchdata
#define FORTgetdata1 getdata1
#define FORTgetdata2 getdata2
#define FORTgetdata2b getdata2b
#define FORTgetdata3 getdata3
#define FORTgetdata2a getdata2a
#define FORTgetsizes getsizes
#define FORTgetsizesa getsizesa
#define FORTgetsizes2 getsizes2
#define FORTgetslicesizes getslicesizes
#define FORTgetslicedata getslicedata
#define FORTgetplot3dq getplot3dq
#define FORTgetplot3dqa getplot3dqa
#define FORTgetsliceparms getsliceparms
#define FORTcolor2rgb color2rgb
#endif

// SUBROUTINE COLOR2RGB(RGB,COLOR)
STDCALL FORTcolor2rgb(int *rgb, char *color, FILE_SIZE colorsize);

STDCALL FORTfcreate_part5sizefile(char *part5file, char *part5sizefile, int *angle_flag, int *error,
                                  FILE_SIZE lenpart5file, FILE_SIZE lenpart5sizefile);

STDCALL FORTgetsliceparms(char *file,int *endian,
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,int *slice3d, int *error,FILE_SIZE lenfile);
STDCALL FORTclosepart(void);
STDCALL FORTclosepatch(void);
STDCALL FORTclosezone(void);
STDCALL FORTgetxyzsize(char *xyzfilename,int *ibp1,int *jbp1,int *kbp1,int *endian,int *error, FILE_SIZE len);
STDCALL FORTgetzonesize(char *zonefilename, int *nzonet, int *nrooms, int *nfires, int *endian, int *error, FILE_SIZE len);
STDCALL FORTgetzonedata(int *nzonet, int *nrooms, int *nfires, 
                        float *zonet, float *zoneqfire, float *zonepr, float *zoneylay,float *zonetl,float *zonetu, int *error);
STDCALL FORTgetxyzdata(int *iblank,int *nx,int *ny,int *nz,int *error);
STDCALL FORTgetpatchsizes1(char *patchfilename,char *patchlonglabel,char *patchshortlabel,char *patchunit,
                           int *endian,int *npatch, int *headersize, int *error,
                           FILE_SIZE len1, FILE_SIZE len2, FILE_SIZE len3, FILE_SIZE len4);
STDCALL FORTgetpatchsizes2(int *version, int *npatch,int *npatchsize, 
                           int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2, int *patchdir,
                           int *headersize, int *framesize);
STDCALL FORTgetpatchdata(int *lunit, int *npatch,int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2,
                         float *patchtimes,float *pqq, int *npqq, int *error);
STDCALL FORTgetdata2a(int *nmax,int *nspr,float *x,float *y,float *z,float *t,
                      float *stime,int *np,int *ns,int *error);
STDCALL FORTgetdata1(int *ipart, int *error);
STDCALL FORTgetdata2(
                     short *xparts, short *yparts, short *zparts,
                     float *t,
                     int *sprinkflag, unsigned char *isprink, float *tspr, int *bframe,int *sframe,int *sprframe,
                     float *ptimes,int *nspr,int *nmax,int *mxframes,int *nframes,
                     int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p, int *partframestep, int *partpointstep, 
                     float *xbox0, float *xbox, float *ybox0, float *ybox, float *zbox0, float *zbox, 
                     float *offset_x, float *offset_y, float *offset_z,
                     int *error,
                     FILE_SIZE lenisprink);
STDCALL FORTgetdata2b(char *partfilename, 
                      float *xpart, float *ypart, float *zpart,  
                      float *tpart,int *bframe,int *sframe,float *ptimes,
                      int *npartpoints,int *npartframes, int *mxframes, int *mxpoints,
                      int *error, FILE_SIZE lenfile);
STDCALL FORTgetdata3(int *nframes,int *npoints,int *error);

STDCALL FORTgetsizesa(char *partfilename, int *npartpoint,int *npartframes,FILE_SIZE lenfile);
STDCALL FORTgetsizes(char *partfilename, int *ibar, int *jbar, int *kbar, 
                     int *nb, int *nv, int *nspr,int *mxframepoints, int *endian, int *showstaticsmoke, int *error, FILE_SIZE filelen);
STDCALL FORTgetsizes2(int *settmin_p, float *tmin_p, int *settmax_p, float *tmax_p,
                      int *nspr, int *frameloadstep, int *partpointstep, int *npartpoints, int *npartframes, int *error);
STDCALL FORTgetslicesizes(char *slicefilename, int *nslicei, int *nslicej, int *nslicek, 
                          int *nsteps,int *sliceframestep, int *endian, int *error, 
                          int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p,
                          int *headersize, int *framesize, int *statfile,
                          FILE_SIZE slicefilelen);
STDCALL FORTgetslicedata(char *slicefilename, char *longlabels, char *shortlabels, 
                         char *units, int *is1,int *is2,int *js1,int *js2, int *ks1,int *ks2,
                         int *idir, float *qslicemin,float *qslicemax,
                         float *qslicedata,float *slicetimes, int *nsteps, int *sliceframestep, int *endian,
                         int *settime_p, int *settmax_p, float *tmin_p, float *tmax_p,
                         FILE_SIZE slicefilelen, FILE_SIZE longlableslen, FILE_SIZE shortlabelslen, FILE_SIZE unitslen);
STDCALL FORTgetplot3dq(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, int *endian, int *isotest, FILE_SIZE filelen);

STDCALL FORTgetplot3dqa(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, FILE_SIZE filelen);

#endif
