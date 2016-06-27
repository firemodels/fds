#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "gd.h"
#include "gd_errors.h"

#ifdef HAVE_LIBVPX
#include "webpimg.h"
#include "gdhelpers.h"

extern void gd_YUV420toRGBA(uint8* Y,
				  uint8* U,
				  uint8* V,
				  gdImagePtr im);

extern void gd_RGBAToYUV420(gdImagePtr im2,
				  uint8* Y,
				  uint8* U,
				  uint8* V);

const char * gdWebpGetVersionString()
{
	return "not defined";
}

BGD_DECLARE(gdImagePtr) gdImageCreateFromWebp (FILE * inFile)
{
	gdImagePtr im;
	gdIOCtx *in = gdNewFileCtx(inFile);
	im = gdImageCreateFromWebpCtx(in);
	in->gd_free(in);

	return im;
}

BGD_DECLARE(gdImagePtr) gdImageCreateFromWebpPtr (int size, void *data)
{
	int	width, height, ret;
 	unsigned char   *Y = NULL;
	unsigned char   *U = NULL;
	unsigned char   *V = NULL;
	gdImagePtr im;

	ret = WebPDecode(data, size, &Y, &U, &V, &width, &height);
	if (ret != webp_success) {
		if (Y) free(Y);
		if (U) free(U);
		if (V) free(V);
		gd_error("WebP decode: fail to decode input data");
		return NULL;
	}
	im = gdImageCreateTrueColor(width, height);
	if (!im) {
		return NULL;
	}
	gd_YUV420toRGBA(Y, U, V, im);
	return im;
}

#define GD_WEBP_ALLOC_STEP (4*1024)

BGD_DECLARE(gdImagePtr) gdImageCreateFromWebpCtx (gdIOCtx * infile)
{
	int	width, height, ret;
	unsigned char   *filedata = NULL;
	unsigned char   *read, *temp;
	unsigned char   *Y = NULL;
	unsigned char   *U = NULL;
	unsigned char   *V = NULL;
	size_t size = 0, n;
	gdImagePtr im;

	do {
		temp = gdRealloc(filedata, size+GD_WEBP_ALLOC_STEP);
		if (temp) {
			filedata = temp;
			read = temp + size;
		} else {
			if (filedata) {
				gdFree(filedata);
			}
			gd_error("WebP decode: realloc failed");
			return NULL;
		}

		n = gdGetBuf(read, GD_WEBP_ALLOC_STEP, infile);
		size += n;
	} while (n>0);

	ret = WebPDecode(filedata, size, &Y, &U, &V, &width, &height);
	gdFree(filedata);
	if (ret != webp_success) {
		if (Y) free(Y);
		if (U) free(U);
		if (V) free(V);
		gd_error("WebP decode: fail to decode input data");
		return NULL;
	}
	im = gdImageCreateTrueColor(width, height);
	gd_YUV420toRGBA(Y, U, V, im);
	return im;
}

BGD_DECLARE(void) gdImageWebpEx (gdImagePtr im, FILE * outFile, int quantization)
{
	gdIOCtx *out = gdNewFileCtx(outFile);
	gdImageWebpCtx(im, out, quantization);
	out->gd_free(out);
}

BGD_DECLARE(void) gdImageWebp (gdImagePtr im, FILE * outFile)
{
	gdIOCtx *out = gdNewFileCtx(outFile);
  	gdImageWebpCtx(im, out, -1);
	out->gd_free(out);
}

BGD_DECLARE(void *) gdImageWebpPtr (gdImagePtr im, int *size)
{
	void *rv;
	gdIOCtx *out = gdNewDynamicCtx(2048, NULL);
	gdImageWebpCtx(im, out, -1);
	rv = gdDPExtractData(out, size);
	out->gd_free(out);

	return rv;
}

BGD_DECLARE(void *) gdImageWebpPtrEx (gdImagePtr im, int *size, int quantization)
{
	void *rv;
	gdIOCtx *out = gdNewDynamicCtx(2048, NULL);
	gdImageWebpCtx(im, out, quantization);
	rv = gdDPExtractData(out, size);
	out->gd_free(out);
	return rv;
}

/*
 * Maps normalized QP (quality) to VP8 QP
 */
int mapQualityToVP8QP(int quality) {
#define MIN_QUALITY 0
#define MAX_QUALITY 100
#define MIN_VP8QP 1
#define MAX_VP8QP 63
	const float scale = MAX_VP8QP - MIN_VP8QP;
	const float vp8qp =
	scale * (MAX_QUALITY - quality) / (MAX_QUALITY - MIN_QUALITY) + MIN_VP8QP;
	if (quality < MIN_QUALITY || quality > MAX_QUALITY) {
		gd_error("Wrong quality value %d.", quality);
		return -1;
	}

	return (int)(vp8qp + 0.5);
}

/* This routine is based in part on code from Dale Lutz (Safe Software Inc.)
 *  and in part on demo code from Chapter 15 of "PNG: The Definitive Guide"
 *  (http://www.cdrom.com/pub/png/pngbook.html).
 */
BGD_DECLARE(void) gdImageWebpCtx (gdImagePtr im, gdIOCtx * outfile, int quantization)
{
	int width = im->sx;
	int height = im->sy;

	int  yuv_width, yuv_height, yuv_nbytes, ret;
	int vp8_quality;
	unsigned char *Y = NULL,
				  *U = NULL,
				  *V = NULL;
	unsigned char *filedata = NULL;

	/* Conversion to Y,U,V buffer */
	yuv_width = (width + 1) >> 1;
	yuv_height = (height + 1) >> 1;
	yuv_nbytes = width * height + 2 * yuv_width * yuv_height;

	if ((Y = (unsigned char *)gdCalloc(yuv_nbytes, sizeof(unsigned char))) == NULL) {
		gd_error("gd-webp error: cannot allocate Y buffer");
		return;
	}
	vp8_quality = mapQualityToVP8QP(quantization);

	U = Y + width * height;
	V = U + yuv_width * yuv_height;
	gd_RGBAToYUV420(im, Y, U, V);

	/* Encode Y,U,V and write data to file */
	ret = WebPEncode(Y, U, V, width, height, width, yuv_width, yuv_height, yuv_width,
					 vp8_quality, &filedata, &yuv_nbytes, NULL);
	gdFree(Y);

	if (ret != webp_success) {
		if (filedata) {
			free(filedata);
		}
		gd_error("gd-webp error: WebP Encoder failed");
		return;
	}

	gdPutBuf (filedata, yuv_nbytes, outfile);
	free(filedata);
}

#endif /* HAVE_LIBVPX */
