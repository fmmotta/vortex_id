#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include "mt64.h"
#include "floodFill.h"
#include "lambdaInit.h"
#include "vortexGen.h"
#include "vortexExtraction.h"

int main(int argc,char **argv){
  const int Width = 100, Height = 100,Pop=10,nVortex=5,nRuns=1000;
  int seed=98755;
  int i,j,err,ngbr,found,nCnect,*label,n,bin,nMax=20;
  int nbList[8],eqList[Pop],**eqClass;
  float Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.;
  float xmin[2]={1.,1.},xmax[2]={9.,9.};
  float *parVortex=NULL,x0[2],dx[2],xf[2],*sField=NULL,*gField;
  float x,y,v0y0 = 0.00,*vCatalog=NULL;
  FILE *dadosgen,*dadosout;
  int hNG=55,hNRc=55,hNa=40,hNb=40,hNN=10;
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  
  x0[0]=0.; xf[0]= 10.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=0.; xf[1]= 10.; dx[1] = (xf[1]-x0[1])/Width;
  
  gField = (float *)malloc(4*Height*Width*sizeof(float));
  if(gField==NULL){
    printf("memory not allocked\n");
    return 1;
  }
  
  sField = (float *)malloc(Height*Width*sizeof(float));
  if(sField==NULL){
    printf("memory not allocked\n");
    return 1;
  }

  label = (int*)malloc(Height*Width*sizeof(int));
  if(label==NULL){
    printf("memory not allocked\n");
    return 2;
  }

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }

  vCatalog = (float*)malloc(4*nMax*sizeof(float));
  if(vCatalog==NULL){
    printf("memory not allocked\n");
    return 3;
  }

  /* histogram preparation - begin */
  hG = gsl_histogram_alloc(hNG); 
  gsl_histogram_set_ranges_uniform(hG,0.,2.5*Gmax);
  hRc = gsl_histogram_alloc(hNRc); 
  gsl_histogram_set_ranges_uniform(hRc,0.,2.5*rmax);
  ha = gsl_histogram_alloc(hNa); 
  gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb = gsl_histogram_alloc(hNb); 
  gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN = gsl_histogram_alloc(hNN);
  gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);
  /* histogram preparation - end*/

  if(argc>1)
    seed = atoi(argv[1]);
  else
    seed = (int) time(NULL);
 
  dadosgen=fopen("data/multiRunGen.txt","w");
  fprintf(dadosgen,"seed: %d\n",seed);
  fprintf(dadosgen,"\ndomain: xi xf dx\n");
  fprintf(dadosgen,"%f %f %f\n",x0[0],xf[0],dx[0]);
  fprintf(dadosgen,"%f %f %f\n",x0[1],xf[1],dx[1]);
  fprintf(dadosgen,"\nvortex params: Binary\n");
  fprintf(dadosgen,"G : %f %f\n",Gmin,Gmax);
  fprintf(dadosgen,"rc: %f %f\n",rmin,rmax);
  fprintf(dadosgen,"a : %f %f\n",xmin[0],xmax[0]);
  fprintf(dadosgen,"b : %f %f\n",xmin[1],xmax[1]);
  fprintf(dadosgen,"\nshear v0/y0=%f\n",v0y0);
  fclose(dadosgen);

  for(n=0;n<nRuns;n+=1){
    //printf("n=%d\n",n);
    err=genLOseenBinaryList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                             nVortex,&parVortex);
    if(err!=0)
      return err;

    for(i=0;i<4*Height*Width;i+=1)
      gField[i]=0.;

    err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0)
      printf("Problems in addSingleOseen\n");

    err=addConstXYShear(x0,dx,Height,Width,v0y0,&gField);
    if(err!=0)
      printf("Problems in addConstXYShear\n");

    err=vortexExtSimple(Height,Width,x0,dx,eqClass,sField,gField,label,
                        0.0,&nCnect,&vCatalog);
    if(err!=0){
      printf("error on vortexExtSimple - %d\n",err);
      return err; 
    }

    gsl_histogram_increment(hN,nCnect);
    for(i=0;i<nCnect;i+=1){
      gsl_histogram_increment(hG,vCatalog[4*i+0]);
      gsl_histogram_increment(hRc,vCatalog[4*i+1]);
      gsl_histogram_increment(ha,vCatalog[4*i+2]);
      gsl_histogram_increment(hb,vCatalog[4*i+3]);
    }

  }

  dadosout=fopen("data/histoG.txt","w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoRc.txt","w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoa.txt","w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histob.txt","w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoN.txt","w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);

  if(sField!=NULL)
    free(sField);
  if(gField!=NULL)
    free(gField);
  if(label!=NULL)
    free(label);
  if(vCatalog!=NULL)
    free(vCatalog);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  free(parVortex); 

  /* histogram free - begin */
  gsl_histogram_free(hG);
  gsl_histogram_free(hRc);
  gsl_histogram_free(ha);
  gsl_histogram_free(hb);
  gsl_histogram_free(hN);
  /* histogram free - end*/

  return 0;
}
