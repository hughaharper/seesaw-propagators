/* imgfftio.c -- routines fimgr and fimgw called by filter_img.f
	filename string must be fixed with '\0' before call.
	fimgw opens for append, so main program should remove old
	file before starting a new run. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void fimgr_(savex, nximg, nyfft, jrowseek, ikillm, ifilnam)
float	*savex;
int	*nximg, *nyfft, *jrowseek, *ikillm;
char	*ifilnam;
{
	int	ntot, jseek, k;
	FILE	*fp = NULL;
	short int	dummy;
	double 	sum;

	sum += 0.0;

	if ((fp = fopen(ifilnam, "r")) == NULL) {
		fprintf(stderr,"fimgr ERROR cannot open r %s\n", ifilnam);
		exit(-1);
	}

	jseek = *jrowseek;
	ntot = (*nximg) * (*nyfft);

	if (jseek && fseek(fp, jseek*(*nximg)*2, 0) ) {
		fprintf(stderr,"fimgr ERROR seeking %d rows in %s\n", jseek,
					ifilnam);
		exit(-1);
	}

	for (k = 0; k < ntot; k++) {
		if ( (fread((char *)&dummy, 2, 1, fp)) != 1) {
			fprintf(stderr,"fimgr ERROR reading from %s\n",
				ifilnam);
			exit(-1);
		}
		savex[k] = (float)dummy;
		sum += (double)dummy;
	}

	fclose(fp);

	if (*ikillm) {
		sum /= ntot;
		for (k = 0; k < ntot; k++) savex[k] -= sum;
	}

	return;
}

void fimgw_(blend, nximg, nyfft2, ofilnam)
float	*blend;
int	*nximg, *nyfft2;
char	*ofilnam;
{
	int	ntot, k, idummy;
	FILE	*fp = NULL;
	short int	dummy;

	if ((fp = fopen(ofilnam, "a")) == NULL) {
		fprintf(stderr,"fimgw ERROR cannot open a %s\n", ofilnam);
		exit(-1);
	}

	ntot = (*nximg) * (*nyfft2);

	for (k = 0; k < ntot; k++) {
		idummy = (int)floor((double)blend[k]+0.5);
		if (abs(idummy) > 32767) {
			fprintf(stderr,"fimgw ERROR data exceeds 32767\n");
			dummy = (idummy < 0) ? -32767 : 32767;
		}
		else {
			dummy = (short)idummy;
		}
		if ( (fwrite((char *)&dummy, 2, 1, fp)) != 1) {
			fprintf(stderr,"fimgw ERROR writing to %s\n",
				ofilnam);
			exit(-1);
		}
	}
	
	fclose(fp);
	return;
}
