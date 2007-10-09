/* $Id: bug00011.c,v 1.1 2007/01/23 23:57:54 pajoye Exp $ */
#include "gd.h"
#include <stdio.h>
#include <stdlib.h>
#include "gdtest.h"

int main()
{
 	gdImagePtr im;
	FILE *fp;

	fp = fopen("emptyfile", "rb");
	if (!fp) {
		printf("failed, cannot open file\n");
	}
	im = gdImageCreateFromPng(fp);
	fclose(fp);

	if (!im) {
		return 0;
	} else {
		return 1;
	}
	//CuAssertTrue(tc, NULL==im);
}
