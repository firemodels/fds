#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "gd.h"
#include "gd_errors.h"

#define TRUE 1
#define FALSE 0

/* Use this for commenting out debug-print statements. */
/* Just use the first '#define' to allow all the prints... */
/*#define GD_SS_DBG(s) (s) */
#define GD_SS_DBG(s)

#ifdef HAVE_LIBPNG
BGD_DECLARE(void) gdImagePngToSink (gdImagePtr im, gdSinkPtr outSink)
{
	gdIOCtx *out = gdNewSSCtx (NULL, outSink);
	gdImagePngCtx (im, out);
	out->gd_free (out);
}

/*
  Function: gdImageCreateFromPngSource

  See <gdImageCreateFromPng> for documentation.  This is obsolete; use
  <gdImageCreateFromPngCtx> instead.
 */
BGD_DECLARE(gdImagePtr) gdImageCreateFromPngSource (gdSourcePtr inSource)
{
	gdIOCtx *in = gdNewSSCtx (inSource, NULL);
	gdImagePtr im;

	im = gdImageCreateFromPngCtx (in);

	in->gd_free (in);

	return im;
}
#else /* no HAVE_LIBPNG */
BGD_DECLARE(void) gdImagePngToSink (gdImagePtr im, gdSinkPtr outSink)
{
	(void)im;
	(void)outSink;
	gd_error("PNG support is not available\n");
}

BGD_DECLARE(gdImagePtr) gdImageCreateFromPngSource (gdSourcePtr inSource)
{
	(void)inSource;
	gd_error("PNG support is not available\n");
	return NULL;
}
#endif /* HAVE_LIBPNG */
