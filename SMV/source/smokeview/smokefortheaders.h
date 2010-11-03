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
#define FORTopenboundary openboundary_
#define FORTfcreate_part5sizefile fcreate_part5sizefile_
#define FORTgetzonesize getzonesize_
#define FORTgetzonedata getzonedata_
#define FORTgetxyzsize getxyzsize_
#define FORTgetxyzdata getxyzdata_
#define FORTgetpatchsizes1 getpatchsizes1_
#define FORTgetpatchsizes2 getpatchsizes2_
#define FORTgetpatchdata getpatchdata_
#define FORTgetdata1 getdata1_
#define FORTgetdata2 getdata2_
#define FORTgetsizes getsizes_
#define FORTgetsizes2 getsizes2_
#define FORTgetsizesa getsizesa_
#define FORTgetslicesizes getslicesizes_
#define FORTgetslicedata getslicedata_
#define FORTgetplot3dq getplot3dq_
#define FORTgetsliceparms getsliceparms_
#define FORTcolor2rgb color2rgb_
#define FORTget_file_unit get_file_unit_
#define FORTclosefortranfile closefortranfile_
#define FORTgetboundaryheader1 getboundaryheader1_
#define FORTgetboundaryheader2 getboundaryheader2_
#else
#define FORTopenboundary openboundary
#define FORTfcreate_part5sizefile fcreate_part5sizefile
#define FORTgetzonesize getzonesize
#define FORTgetzonedata getzonedata
#define FORTgetxyzdata getxyzdata
#define FORTgetpatchsizes1 getpatchsizes1
#define FORTgetpatchsizes2 getpatchsizes2
#define FORTgetpatchdata getpatchdata
#define FORTgetdata1 getdata1
#define FORTgetdata2 getdata2
#define FORTgetsizes getsizes
#define FORTgetsizesa getsizesa
#define FORTgetsizes2 getsizes2
#define FORTgetslicesizes getslicesizes
#define FORTgetslicedata getslicedata
#define FORTgetplot3dq getplot3dq
#define FORTgetsliceparms getsliceparms
#define FORTcolor2rgb color2rgb
#define FORTget_file_unit get_file_unit
#define FORTclosefortranfile closefortranfile 
#define FORTgetboundaryheader1 getboundaryheader1
#define FORTgetboundaryheader2 getboundaryheader2
#endif

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
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,int *slice3d, int *error,FILE_SIZE lenfile);
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
