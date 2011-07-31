// $Date$ 
// $Revision$
// $Author$

#ifdef IN_UPDATE
#define UPEXTERN
#else
#define UPEXTERN extern CCC
#endif

UPEXTERN int compare_float( const void *arg1, const void *arg2 );
UPEXTERN void update_framenumber(int changetime);
UPEXTERN void reset_itimes0(void);
UPEXTERN void updateShow(void);
UPEXTERN void synctimes(void);
UPEXTERN void updatetimes(void);
UPEXTERN int getplotstate(int choice);
UPEXTERN int getindex(float key, const float *list, int nlist);
UPEXTERN int isearch(float *list, int nlist, float key, int guess);




