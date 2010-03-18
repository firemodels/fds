// $Date$ 
// $Revision$
// $Author$

#ifdef pp_WIN_INTEL
#pragma warning (disable:2259)
#pragma warning (disable:1572)
#endif

/* --------------------------  flowlabels ------------------------------------ */

#define NHIST_BUCKETS 100000
typedef struct {
  int buckets[NHIST_BUCKETS];
  float valmin, valmax;
  int ntotal;
} histogramdata;

//************************** headers ****************************************

void init_histogram(histogramdata *histogram);
void copy_data2histogram(float *vals, int nvals, histogramdata *histgram);
void update_histogram(float *vals, int nvals, histogramdata *histogram);
void merge_histogram(histogramdata *histogram1, histogramdata *histogram2);
float get_histogram_value(histogramdata *histogram, float cdf);
void check_histogram(void);


