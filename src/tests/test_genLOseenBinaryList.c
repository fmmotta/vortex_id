#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "vortexGen.h"
#include "mt64.h"

int main(int argc,char **argv){
  const int numVortex=5;
  int i,err,seed=98755;
  double Gmin=1,Gmax=2,rmin=1,rmax=2;
  double xmin[2]={-4.,-4.},xmax[2]={4.,4.};
  double *parVortex=NULL;

  if(argc>1)
    seed = atoi(argv[1]);

  err=genLOseenBinaryList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                           numVortex,&parVortex);

  for(i=0;i<numVortex;i+=1)
    printf("G=%lf rc = %lf a=%lf b=%lf\n",parVortex[4*i+0],
                                        parVortex[4*i+1],
                                        parVortex[4*i+2],
                                        parVortex[4*i+3]);

  free(parVortex); parVortex=NULL;
  
  return 0;
}