#ifndef LIGHTSMOKE_H_DEFINED
#define LIGHTSMOKE_H_DEFINED
//***********************
//************* #definess
//***********************
#ifdef INLIGHTSMOKE
#define LIGHTEXTERN
#else
#define LIGHTEXTERN extern
#endif

//***********************
//************* structures
//***********************

/* --------------------------  radiancedata ------------------------------------ */

typedef struct {
  int *ijkbar;
  unsigned char *radiance, *opacity;
  float *xyzbar0, *xyzbar;
  float *dxyz;
} radiancedata;



//***********************
//************* headers
//***********************

void setup_radiancemap(radiancedata *radianceinfo, int ijkbar[3], float xyzbar0[3], float xyzbar[3], float dxyz[3],
                               unsigned char *radiance, unsigned char *opacity);
void build_radiancemap2(radiancedata *radianceinfo);
void build_radiancemap(radiancedata *radianceinfo);
#ifdef pp_KDTEST
void test_kd(void);
#endif
#endif

