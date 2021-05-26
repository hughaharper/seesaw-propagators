/************************************************************************
* readgrd routine to read an srtm mask and return 1-land or 0-water     *
*************************************************************************/
/************************************************************************
* Creator: David T. Sandwell    Scripps Institution of Oceanography     *
* Date   : 02/12/07                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*   05/16/08:  Karen alters function to be "void" and return idry as an *
*   argument so that the function works with her fortran programs.      *
************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

unsigned char maskset[8] = { 128, 64, 32, 16, 8, 4, 2, 1};

void readsrtmmask_(rln,rlt,idry)

double *rln;		/* longitude -180 to 180 */
double *rlt;		/* latitude -90 to 90  */
int *idry;

{
   char fname[120]="/cryosat/mgg/fftfault0/img/srtm.bitmask";
   FILE *fp;
   static unsigned char *a;
   static int ifirst = 1;
   int i, j, kk, nmask = 116640000;

/* if this is the first call to this routine then malloc the memory 
   and read the bitmask */

   if(ifirst == 1) {
      if(( a = (unsigned char *) malloc(nmask)) == NULL){
           printf("sorry, couldn't allocate memory for mask\n");
           exit(-1);
      }
   
      if ( (fp = fopen(fname, "r")) == NULL) {
           fprintf(stderr,"read_srtm_mask:  ERROR.  Cannot open w %s\n", fname);
           exit(-1);
      }
      if ((fread(a, 1, nmask, fp)) != nmask) {
           fprintf(stderr,"read_srtm_mask:  ERROR during read of %s\n", fname);
           exit(-1);
      }
      ifirst = 0;
   }

/* compute the index in the mask */

   j = (int)floor((*rln + 180.)*120.);
   i = (int)floor((90. - *rlt )*120.);
   kk = 120*360*i + j;
   *idry = 0;
   if(a[kk/8] & maskset[kk%8]) *idry = 1;

/*   printf("lon is %f, lat is %f, idry is %i\n",*rln,*rlt,*idry);
     printf("%i\n",*idry);*/

   return;
}
