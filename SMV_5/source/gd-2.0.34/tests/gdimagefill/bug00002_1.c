#include <gd.h>
#include <stdio.h>
#include "gdtest.h"

int main()
{
	gdImagePtr im;
	int error = 0;
	FILE *fp;

	fputs("flag 0\n", stdout);
	im = gdImageCreateTrueColor(100, 100);

	fputs("flag 1\n", stdout);
	gdImageFill(im, 0, 0, 0xffffff);
	fputs("flag 2\n", stdout);
	gdImageFill(im, 0, 0, 0xffffff);
	fputs("flag 3\n", stdout);

	fp = fopen("a.png", "wb");
 	/* Write img to stdout */
 	gdImagePng(im,fp);
	fclose(fp);

	if (!gdAssertImageEqualsToFile("bug00002_1_exp.png", im)) {
		error = 1;
	}

 	/* Destroy it */
 	gdImageDestroy(im);
	return error;
}

