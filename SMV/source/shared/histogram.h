#ifndef HISTOGRAM_H_DEFINED
#define HISTOGRAM_H_DEFINED
#ifdef pp_INTEL
#pragma warning (disable:2259)
#pragma warning (disable:1572)
#endif

/* --------------------------  flowlabels ------------------------------------ */

#define NHIST_BUCKETS 100000
typedef struct {
  int *buckets, *buckets_2d;
  float *rvals;
  int nbuckets, ndim, defined;
  int nx, ny, ntotal;
  float valmin, valmax;
  float valxmin, valxmax;
  float valymin, valymax;
  int complete;
} histogramdata;

//************************** headers ****************************************

void complete_histogram(histogramdata *histogram);
void copy_data2histogram(float *vals, int nvals, histogramdata *histogram);
void copy_uvdata2histogram(float *uvals, float *vvals, int nvals, histogramdata *histogram);
void free_histogram(histogramdata *histogram);
void free_histogram2d(histogramdata *histogram);
float get_histogram_value(histogramdata *histogram, float cdf);
void init_histogram(histogramdata *histogram, int nbuckets);
void init_histogram2d(histogramdata *histogram, int nx, int ny);
void merge_histogram(histogramdata *histogramto, histogramdata *histogramfrom);
void merge_uvhistogram(histogramdata *histogramto, histogramdata *histogramfrom);
void reset_histogram(histogramdata *histogram);
void update_uvhistogram(float *uvals, float *vvals, int nvals, histogramdata *histogramto);
void update_histogram(float *vals, int nvals, histogramdata *histogram);
#if pp_CHECK
void check_histogram(void);
#endif
#endif
