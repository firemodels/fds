// $Date$ 
// $Revision$
// $Author$

#ifdef pp_WIN_INTEL
#pragma warning (disable:2259)
#pragma warning (disable:1572)
#endif

/* --------------------------  flowlabels ------------------------------------ */

#define NBUCKETS 100000
typedef struct {
  int buckets[NBUCKETS];
  float valmin, valmax;
  int ntotal;
} histogramdata;

//************************** headers ****************************************

void init_histogram(float *vals, int nvals, histogramdata *histogram);
void update_histogram(float *vals, int nvals, histogramdata *histogram);
void merge_histogram(histogramdata *histogram1, histogramdata *histogram2);
float get_histogram_value(histogramdata *histogram, float cdf);
void check_histogram(void);


