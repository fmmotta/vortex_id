#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

int main(int argc,char **argv){
  const int Nbins=100;
  int i,k,err;
  int numG[Nbins+1],numRc[Nbins+1];
  double y0=0.,y1=1.0,dy=0.;
  double G,rc,x,y,M[2][2],rg,L[2],v[2][2],sinRef,S[2][2],omega,strain;
  double avgG[Nbins+1],avgRc[Nbins+1];
  FILE *dados;
  
  dy = (y1-y0)/Nbins;

  dados = fopen(argv[1],"r");
  if(dados==NULL)
    printf("Problems\n");

  for(k=0;k<=Nbins;k+=1){
    numG[k] = 0 ; numRc[k] = 0 ;
    avgG[k] = 0.; avgRc[k] = 0.;
  }

  while( true ){
    err=fscanf(dados,"%lf %lf %lf %lf ", &G,&rc,&x,&y);
    if(err!=4)
      break;
    
    k = (int)( (y-y0)/dy );
    numG[k] += 1; numRc[k] += 1;
    avgG[k] += G; avgRc[k] += rc;
  }

  fclose(dados);

  fopen("results.dat","w");
  for(k=0;k<Nbins;k+=1){
    if(numG[k]!=0)
      avgG[k] /= numG[k];
    if(numRc[k]!=0)
      avgRc[k] /= numRc[k];

    fprintf(dados,"%f %f %f\n",y0+k*dy,avgG[k],avgRc[k]);
  }
  fclose(dados);


  return 0;
}