#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

int main(int argc,char **argv){
  const int Nbins=100;
  int i,k,err;
  int numGp[Nbins+1],numRcp[Nbins+1],numGn[Nbins+1],numRcn[Nbins+1];
  double y0=0.,y1=1.0,dy=0.;
  double G,rc,x,y;
  double avgGp[Nbins+1],avgRcp[Nbins+1],avgGn[Nbins+1],avgRcn[Nbins+1];
  double avgWp[Nbins+1],avgWn[Nbins+1];
  double avgG2p[Nbins+1],avgRc2p[Nbins+1],avgG2n[Nbins+1],avgRc2n[Nbins+1];
  double avgW2p[Nbins+1],avgW2n[Nbins+1];
  FILE *dados;
  
  dy = (y1-y0)/Nbins;

  dados = fopen(argv[1],"r");
  if(dados==NULL)
    printf("Problems\n");

  for(k=0;k<=Nbins;k+=1){
    numGp[k] = 0 ; numRcp[k] = 0 ;
    avgGp[k] = 0.; avgRcp[k] = 0.;
    avgGn[k] = 0.; avgRcn[k] = 0.;
    avgWn[k] = 0.; avgWp[k] = 0.;
    avgG2p[k] = 0.; avgRc2p[k] = 0.;
    avgG2n[k] = 0.; avgRc2n[k] = 0.;
    avgW2n[k] = 0.; avgW2p[k] = 0.;
  }

  while( true ){
    err=fscanf(dados,"%lf %lf %lf %lf ", &G,&rc,&x,&y);
    if(err!=4)
      break;
    
    k = (int)( (y-y0)/dy );
    if(G>=0){
      numGp[k] += 1; numRcp[k] += 1;
      avgGp[k] += G; avgRcp[k] += rc;
      avgWp[k] += G/(rc*rc);

      avgG2p[k] += G*G; avgRc2p[k] += rc*rc;
      avgW2p[k] += (G/(rc*rc))*(G/(rc*rc));
    }
    else{
      numGn[k] += 1; numRcn[k] += 1;
      avgGn[k] += G; avgRcn[k] += rc;
      avgWn[k] += G/(rc*rc);

      avgG2n[k] += G*G; avgRc2n[k] += rc*rc;
      avgW2n[k] += (G/(rc*rc))*(G/(rc*rc));
    }
  }

  fclose(dados);

  dados = fopen(argv[2],"w");
  if(dados==NULL)
  	printf("Shit!\n");
  for(k=0;k<Nbins;k+=1){
    if(numGp[k]!=0)
      avgGp[k] /= numGp[k];
    if(numGp[k]!=0)
      avgWp[k] /= numGp[k];
    if(numRcp[k]!=0)
      avgRcp[k] /= numRcp[k];

    if(numGp[k]!=0)
      avgG2p[k] /= numGp[k];
    if(numGp[k]!=0)
      avgW2p[k] /= numGp[k];
    if(numRcp[k]!=0)
      avgRc2p[k] /= numRcp[k];

    if(numGn[k]!=0)
      avgGn[k] /= numGn[k];
    if(numGn[k]!=0)
      avgWn[k] /= numGn[k];
    if(numRcn[k]!=0)
      avgRcn[k] /= numRcn[k];

    if(numGn[k]!=0)
      avgG2n[k] /= numGn[k];
    if(numGn[k]!=0)
      avgW2n[k] /= numGn[k];
    if(numRcn[k]!=0)
      avgRc2n[k] /= numRcn[k];

    fprintf(dados,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n",y0+k*dy
                                          ,avgGp[k]
                                          ,sqrt(avgG2p[k]-avgGp[k]*avgGp[k])
                                          ,avgRcp[k]
                                          ,sqrt(avgRc2p[k]-avgRcp[k]*avgRcp[k])
                                          ,avgWp[k]
                                          ,sqrt(avgW2p[k]-avgWp[k]*avgWp[k])
                                          ,avgGn[k]
                                          ,sqrt(avgG2n[k]-avgGn[k]*avgGn[k])
                                          ,avgRcn[k]
                                          ,sqrt(avgRc2n[k]-avgRcn[k]*avgRcn[k])
                                          ,avgWn[k]
                                          ,sqrt(avgW2n[k]-avgWn[k]*avgWn[k]));
  }
  fclose(dados);


  return 0;
}