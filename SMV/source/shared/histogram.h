#ifndef HISTOGRAM_H_DEFINED
#define HISTOGRAM_H_DEFINED
#ifdef pp_INTEL
#pragma warning (disable:2259)
#pragma warning (disable:1572)
#endif

#define HIST_USE_BOUNDS 0
#define HIST_COMPUTE_BOUNDS 1
#define HIST_OK 0
#define HIST_OLD 1
#define HIST_ERR -1


/* --------------------------  flowlabels ------------------------------------ */

#define NHIST_BUCKETS 100000
typedef struct {
  int *buckets, *buckets_2d;
  int nbuckets, ndim, defined;
  int nx, ny, ntotal;
  float valmin, valmax;
  float valxmin, valxmax;
  float valymin, valymax;
  float valmean, valstdev;
  int complete;
} histogramdata;

//************************** headers ****************************************

void complete_histogram(histogramdata *histogram);
void copy_buckets2histogram(int *buckets, int nbuckets, float valmin, float valmax, histogramdata *histogram);
void copy_data2histogram(float *vals, int nvals, histogramdata *histogram);
void copy_polardata2histogram(float *speed, float *angle, int nvals, float rmin, float rmax, histogramdata *histogram);
void copy_uvdata2histogram(float *uvals, float *vvals, int nvals, float rmin, float rmax, histogramdata *histogram);

void free_histogram(histogramdata *histogram);
void free_histogram2d(histogramdata *histogram);
void get_2dminmax(float *uvals, float *vvals, int nvals, float *rmin, float *rmax, int flag);
void get_polarminmax(float *speed, int nvals, float *rmin, float *rmax, int flag);
float get_histogram_value(histogramdata *histogram, float cdf);
void get_histogram_statistics(histogramdata *histogram);
void init_histogram(histogramdata *histogram, int nbuckets);
void init_histogram2d(histogramdata *histogram, int nx, int ny);
void merge_histogram(histogramdata *histogramto, histogramdata *histogramfrom);
void reset_histogram(histogramdata *histogram);
void reset_histogram2d(histogramdata *histogram);
void update_histogram(float *vals, int nvals, histogramdata *histogram);

#endif
