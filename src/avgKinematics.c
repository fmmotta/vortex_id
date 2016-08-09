#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

int main(int argc,char **argv){
  const int Nbins=100;
  int k,err;
  int numGp[Nbins+1],numRcp[Nbins+1],numGn[Nbins+1],numRcn[Nbins+1],count=0;
  double y0=0.,y1=1.0,dy=0.;
  double G,rc,x,y,u,v;
  double avgGp[Nbins+1],avgRcp[Nbins+1],avgGn[Nbins+1],avgRcn[Nbins+1];
  double avgWp[Nbins+1],avgWn[Nbins+1];
  double avgUp[Nbins+1],avgUn[Nbins+1],avgVp[Nbins+1],avgVn[Nbins+1];
  FILE *dados;
  
  dy = (y1-y0)/Nbins;

  dados = fopen(argv[1],"r");
  if(dados==NULL)
    printf("Problems\n");

  for(k=0;k<=Nbins;k+=1){
    numGp[k] = 0 ; numRcp[k] = 0 ;
    avgGp[k] = 0.; avgRcp[k] = 0.;
    numGn[k] = 0 ; numRcn[k] = 0 ;
    avgGn[k] = 0.; avgRcn[k] = 0.;
    avgWp[k] = 0.; avgWn[k] = 0.;
    avgUp[k] = 0.; avgUn[k] = 0.;
    avgVp[k] = 0.; avgVn[k] = 0.;
  }

  count=0;
  while( true ){
    err=fscanf(dados,"%lf %lf %lf %lf %lf %lf", &G,&rc,&x,&y,&u,&v);
    if(err!=6)
      break;
    
    k = (int)( (y-y0)/dy );
    if(G>=0){
      numGp[k] += 1; numRcp[k] += 1;
      avgGp[k] += G; avgRcp[k] += rc;
      avgWp[k] += G/(rc*rc);
      avgUp[k] += u;
      avgVp[k] += v;
    }
    else{
      numGn[k] += 1; numRcn[k] += 1;
      avgGn[k] += G; avgRcn[k] += rc;
      avgWn[k] += G/(rc*rc);
      avgUn[k] += u;
      avgVn[k] += v;
    }
  }

  fclose(dados);

  dados=fopen(argv[2],"w");
  for(k=0;k<Nbins;k+=1){
    if(numGp[k]!=0){
      avgGp[k] /= numGp[k];
      avgWp[k] /= numGp[k];
      avgRcp[k] /= numGp[k];
      avgUp[k] /= numGp[k];
      avgVp[k] /= numGp[k];
    }

    if(numGn[k]!=0){
      avgGn[k] /= numGn[k];
      avgWn[k] /= numGn[k];
      avgRcn[k] /= numGn[k];
      avgUn[k] /= numGn[k];
      avgVn[k] /= numGn[k];
    }

    fprintf(dados,"%f %f %f %f %f %f %f %d %d\n",y0+k*dy,avgGp[k],avgRcp[k]
                                                ,avgGn[k],avgRcn[k]
                                                ,avgWp[k],avgWn[k]
                                                ,numGp[k],numGn[k]);
  }
  fclose(dados);

  if(argc>=4){
    FILE *dadosKin=fopen(argv[3],"w");
    for(k=0;k<Nbins;k+=1){
      fprintf(dadosKin,"%lf %lf %lf %f %f\n",y0+k*dy,avgUn[k],avgUp[k]
                                                    ,avgVn[k],avgVp[k]);
    }
    fclose(dadosKin);
  }

  return 0;
}