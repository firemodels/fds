/* $Id: bug00020.c,v 1.1 2007/01/23 23:57:53 pajoye Exp $ */
#include "gd.h"
#include "gdtest.h"

#define exp_img "bug00020_exp.png"
#define width 50

int main()
{
 	gdImagePtr im, im2;
 	int error = 0;

	im = gdImageCreateTrueColor(width, width);
	gdImageFilledRectangle(im, 0,0, width, width, 0xFF0000);
	gdImageColorTransparent(im, 0xFF0000);
	gdImageFilledEllipse(im, width/2, width/2, width - 20, width - 10,
				0x50FFFFFF);

	im2 = gdImageCreateTrueColor(width, width);

	gdImageCopyRotated(im2, im, width / 2, width / 2, 0,0, width, width, 60);

	if (!gdAssertImageEqualsToFile(exp_img, im2)) {
		error = 1;
	}

	gdImageDestroy(im2);
 	gdImageDestroy(im);
 	return error;
}
