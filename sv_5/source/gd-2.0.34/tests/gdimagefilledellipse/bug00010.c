/* $Id: bug00010.c,v 1.1 2007/01/23 23:57:53 pajoye Exp $ */
#include "gd.h"
#include "gdtest.h"

#define exp_img "bug00010_exp.png"

int main()
{
 	gdImagePtr im;
 	int error = 0;

	im = gdImageCreateTrueColor(100,100);
	gdImageFilledEllipse(im, 50,50, 70, 90, 0x50FFFFFF);
	if (!gdAssertImageEqualsToFile(exp_img, im)) {
		error = 1;
	}

 	gdImageDestroy(im);
 	return error;
}
