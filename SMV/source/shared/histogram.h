#ifndef HISTOGRAM_H_DEFINED
#define HISTOGRAM_H_DEFINED
#ifdef pp_INTEL
#pragma warning (disable:2259)
#pragma warning (disable:1572)
#endif

/* --------------------------  flowlabels ------------------------------------ */

#define NHIST_BUCKETS 100000
typedef struct {
  int buckets[NHIST_BUCKETS], *buckets2;
  float *rvals;
  int nbuckets, ndim, defined;
  int nx, ny, ntotal;
  float valmin, valmax;
  float valxmin, valxmax;
  float valymin, valymax;
  int complete;
} histogramdata;

//************************** headers ****************************************

void init_histogram(histogramdata *histogram);
void init_histogram2d(histogramdata *histogram, int nx, int ny);
void copy_data2histogram(float *vals, int nvals, histogramdata *histgram);
void copy_uvdata2histogram(float *uvals, float *vvals, int nvals, histogramdata *histgram);
void update_histogram(float *vals, int nvals, histogramdata *histogram);
void merge_histogram(histogramdata *histogram1, histogramdata *histogram2);
float get_histogram_value(histogramdata *histogram, float cdf);
void check_histogram(void);
void complete_histogram(histogramdata *histgram);
#endif
