// $Date$ 
// $Revision$
// $Author$
#ifndef UPDATE_H_DEFINED
#define UPDATE_H_DEFINED
#ifdef IN_UPDATE
#define UPEXTERN
#else
#define UPEXTERN extern CCC
#endif

UPEXTERN void updateclip(int slicedir);
UPEXTERN void updateclipbounds(int set_i0, int *i0, int set_i1, int *i1, int maxi);
UPEXTERN int compare_float( const void *arg1, const void *arg2 );
UPEXTERN void update_framenumber(int changetime);
UPEXTERN void reset_itimes0(void);
UPEXTERN void updateShow(void);
UPEXTERN void synctimes(void);
UPEXTERN void updatetimes(void);
UPEXTERN int getplotstate(int choice);
UPEXTERN int getindex(float key, const float *list, int nlist);
UPEXTERN int isearch(float *list, int nlist, float key, int guess);
#endif



