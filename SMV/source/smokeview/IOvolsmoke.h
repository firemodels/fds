#ifndef IOVOLSMOKE_H_DEFINED
#define IOVOLSMOKE_H_DEFINED
#define GPU_VOLframemax 1.5
#define MLEFT 0
#define MFRONT 1
#define MDOWN 2
#define MRIGHT 3
#define MBACK 4
#define MUP 5
#define MEPS 0.1

#ifdef pp_GPU
void init_volsmoke_supertexture(supermeshdata *smesh);
void init_volsmoke_texture(meshdata *meshi);
void update_volsmoke_texture(meshdata *meshi, float *smokedata, float *firedata);
#endif
void init_supermesh(void);
void unload_volsmoke_frame_allmeshes(int framenum);
void compute_all_smokecolors(void);
void drawsmoke3dGPUVOL(void);
void drawsmoke3dVOL(void);
void drawsmoke3dVOLdebug(void);
void get_cum_smokecolor(float *cum_smokecolor, float *xyzvert, float dstep, meshdata *meshi, int iwall);
void get_pt_smokecolor(float *smoke_tran, float **smoke_color, float dstep, float xyz[3], meshdata *meshi, int *inobst, char *blank);
void get_volsmoke_all_times(volrenderdata *vr);
float get_volsmoke_frame_time(volrenderdata *vr, int framenum);
int get_volsmoke_nframes(volrenderdata *vr);
void init_volrender(void);
void init_volrender_surface(int firstcall);
void read_volsmoke_allframes(volrenderdata *vr);
void read_volsmoke_allframes_allmeshes(void);
void free_volsmoke_frame(volrenderdata *vr, int framenum);
void read_volsmoke_frame(volrenderdata *vr, int framenum, int *first);
void read_volsmoke_frame_allmeshes(int framenum, supermeshdata *smesh);
void unload_volsmoke_allframes(volrenderdata *vr);
void *read_volsmoke_allframes_allmeshes2(void *arg);
#endif

