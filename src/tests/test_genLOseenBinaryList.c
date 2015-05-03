#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../vortexGen.h"
#include "../mt64.h"

int main(int argc,char **argv){
  const int numVortex=5;
  int i,err,seed=98755;
  float Gmin=1,Gmax=2,rmin=1,rmax=2;
  float xmin[]={-4.,-4.},xmax[]={4.,4.};
  float *parVortex=NULL;

  if(argc>1)
    seed = atoi(argv[1]);

  err=genLOseenBinaryList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                           numVortex,&parVortex);

  for(i=0;i<numVortex;i+=1)
    printf("G=%lf rc = %lf a=%lf b=%lf\n",parVortex[4*i+0],
                                        parVortex[4*i+1],
                                        parVortex[4*i+2],
                                        parVortex[4*i+3]);

  return 0;
}