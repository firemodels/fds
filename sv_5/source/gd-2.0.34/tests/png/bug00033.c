/* $Id: bug00033.c,v 1.1 2007/01/23 23:57:54 pajoye Exp $ */
#include "gd.h"
#include <stdio.h>
#include <stdlib.h>
#include "gdtest.h"

int main()
{
 	gdImagePtr im;
	FILE *fp;

	fp = fopen("bug00033.png", "rb");
	if (!fp) {
		printf("failed, cannot open file <bug00033.png>\n");
		return 1;
	}

	im = gdImageCreateFromPng(fp);
	fclose(fp);

	if (im) {
		gdImageDestroy(im);
		return 1;
	} else {
		return 0;
	}

}
