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

void CompleteHistogram(histogramdata *histogram);
void CopyBuckets2Histogram(int *buckets, int nbuckets, float valmin, float valmax, histogramdata *histogram);
void CopyU2Histogram(float *vals, int nvals, histogramdata *histogram);
void CopyPolar2Histogram(float *speed, float *angle, int nvals, float rmin, float rmax, histogramdata *histogram);
void CopyUV2Histogram(float *uvals, float *vvals, int nvals, float rmin, float rmax, histogramdata *histogram);

void FreeHistogram(histogramdata *histogram);
void FreeHistogram2d(histogramdata *histogram);
void Get2DMinMax(float *uvals, float *vvals, int nvals, float *rmin, float *rmax, int flag);
void GetPolarMinMax(float *speed, int nvals, float *rmin, float *rmax, int flag);
float GetHistogramVal(histogramdata *histogram, float cdf);
void GetHistogramStats(histogramdata *histogram);
void InitHistogram(histogramdata *histogram, int nbuckets);
void InitHistogram2D(histogramdata *histogram, int nx, int ny);
void MergeHistogram(histogramdata *histogramto, histogramdata *histogramfrom);
void ResetHistogram(histogramdata *histogram);
void ResetHistogram2d(histogramdata *histogram);
void UpdateHistogram(float *vals, int nvals, histogramdata *histogram);

#endif
