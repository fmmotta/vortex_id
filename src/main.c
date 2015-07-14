#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_histogram.h>
#include "ini.h"
#include "mt64.h"
#include "vortexGen.h"
#include "floodFill.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "inputManager.h"

#define DEBUG_MODE false

int fprintfRunParamSigned(FILE *dadosgen,long long int seed,float x0[],
                         float xf[],float dx[],float Gmin, float Gmax,
                         float rmin, float rmax, float xmin[], float xmax[], 
                         float v0y0);

int histoIncVortex(int nVortex, float *parVortex,
                   gsl_histogram *iG, gsl_histogram *iRc,
                   gsl_histogram *ia, gsl_histogram *ib);

int fprintVortex(FILE *dadosout, int run,int nVortex, float *vCatalog);

int fprintsField(FILE *dadosout,float *x0,float *dx,
                 int Width, int Height, float *sField);

int fprintLabels(FILE *dadosout,float *x0,float *dx,
                 int Width, int Height, int *label);

int main(int argc,char **argv){
  int Width = 100, Height = 100,nVortex=5,nFixVortex=5,nRuns=1000;
  int runType=0;
  int numG=3,numRc=3;
  int *label=NULL,**eqClass=NULL;
  long long int seed=98755;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10;
  int i,j,err,nCnect,rCnect=0,n,it,nMax=500,pass=0;
  float Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  float xmin[2]={-9.,-9.},xmax[2]={9.,9.};
  float *parVortex=NULL,*Glist,*Rclist;
  float x0[2],dx[2],xf[2],*sField=NULL,*gField=NULL;
  float x,y,v0y0 = 0.00,*vCatalog=NULL,*rCatalog=NULL,*majorVortex=NULL;
  float hGmin=0.,hGmax=0.,hRcMin=0.,hRcMax=0.;
  char genFile[300+1],folder[100+1],tag[100+1],filename[400+1];
  FILE *dadosgen,*dadosout,*dadosVin,*dadosVout;
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  gsl_histogram *iG,*iRc,*ia,*ib;
  configVar cfg;
  
  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly the configuration file\n");
    return -1;
  }

  /* Loading Configuration -- I need something more concise */

  err=initConfig(&cfg);

  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
  
  if(DEBUG_MODE==true){
    err=printConfig(&cfg);
    if(err!=0)
      return err;
  }
  
  if(cfg.dim!=2){
    printf("Dimension is not 2 - Can't Follow\n");
    return 2;
  }

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }

  seed    = cfg.seed;
  Width   = cfg.Width;
  Height  = cfg.Height;
  nRuns   = cfg.nRuns;
  runType = cfg.runType;
  nFixVortex = cfg.nVortex;

  numG    = cfg.numG;
  numRc   = cfg.numRc;

  Glist   = (float*)malloc(numG*sizeof(float));
  if(Glist==NULL){printf("Can't allocate Glist\n"); return 3;}
  for(i=0;i<numG;i+=1)
    Glist[i] = cfg.Glist[i];

  Rclist   = (float*)malloc(numRc*sizeof(float));
  if(Rclist==NULL){printf("Can't allocate Glist\n"); return 3;}
  for(i=0;i<numRc;i+=1)
    Rclist[i]=cfg.Rclist[i];

  x0[0]   = cfg.x0[0]; x0[1]   = cfg.x0[1]; 
  xmin[0] = x0[0]+1; xmin[1] = x0[1]+1;
  
  xf[0]   = cfg.xf[0]; xf[1]   = cfg.xf[1]; 
  xmax[0] = xf[0]-1; xmax[1] = xf[1]-1;
  
  dx[0] = (xf[0]-x0[0])/Height;
  dx[1] = (xf[1]-x0[1])/Width;
  
  Gmin    = cfg.Gmin;
  Gmax    = cfg.Gmax;
  rmin    = cfg.RcMin;
  rmax    = cfg.RcMax;
  v0y0    = cfg.v0y0;

  if(runType==0)
    threshold = cfg.swThresh;
  else if (runType==2)
    threshold = cfg.sndSwThresh;
  else
    threshold = -1.;

  hNG     = cfg.hNG;
  hNRc    = cfg.hNRc;
  hNa     = cfg.hNa;
  hNb     = cfg.hNb;
  hNN     = cfg.hNN;
  hGmin   = cfg.hGmin;
  hGmax   = cfg.hGmax;
  hRcMin  = cfg.hRcMin;
  hRcMax  = cfg.hRcMax;
  
  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */
  /* Memory Allocation */

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
  gsl_histogram_set_ranges_uniform(hG,hGmin,hGmax);
  hRc = gsl_histogram_alloc(hNRc); 
  gsl_histogram_set_ranges_uniform(hRc,hRcMin,hRcMax);
  ha = gsl_histogram_alloc(hNa); 
  gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb = gsl_histogram_alloc(hNb); 
  gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN = gsl_histogram_alloc(hNN);
  gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);

  iG = gsl_histogram_alloc(hNG); 
  gsl_histogram_set_ranges_uniform(iG,hGmin,hGmax);
  iRc = gsl_histogram_alloc(hNRc); 
  gsl_histogram_set_ranges_uniform(iRc,hRcMin,hRcMax);
  ia = gsl_histogram_alloc(hNa); 
  gsl_histogram_set_ranges_uniform(ia,xmin[0],xmax[0]);
  ib = gsl_histogram_alloc(hNb); 
  gsl_histogram_set_ranges_uniform(ib,xmin[1],xmax[1]);
  /* histogram preparation - end*/

  sprintf(genFile,"%s/genfile-%s.dat",folder,tag);
  dadosgen=fopen(genFile,"w");
  err=fprintfRunParamSigned(dadosgen,seed,x0,xf,dx,Gmin,Gmax,rmin,
                            rmax,xmin,xmax,v0y0);
  fclose(dadosgen);

  sprintf(filename,"%s/inputVortexes.txt",folder);
  dadosVin = fopen(filename,"w");
  sprintf(filename,"%s/outputVortexes.txt",folder);
  dadosVout = fopen(filename,"w");

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  for(n=0;n<nRuns;n+=1){
    if(n%1000 == 0){
      printf("%d runs have passed\n",n);
      fflush(dadosVin);
      fflush(dadosVout);
    }

    nVortex = nFixVortex;
    err=genLOseenSignUniformList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                                 nVortex,&parVortex);
    if(err<0)
      return err;
    else if((err>0) && (err<nVortex))
      nVortex = err;

    for(i=0;i<4*Height*Width;i+=1)
      gField[i]=0.;

    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;

    err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0)
      printf("Problems in addSingleOseen\n");

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0)
      printf("Problems in gradUtoLamb\n");
  
    err = floodFill(sField,Width,Height,eqClass,label);
    if(err!=0)
      printf("Problems in floodFill\n");

    err = renameLabels(Height,Width,label);
    if(err>0)
      nCnect=err;
    else
      printf("problems with renameLabels - %d\n",err);

    err=vortexExtraction(Height,Width,nCnect,x0,dx,sField,
                         gField,label,&vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }

    err=histoIncVortex(nVortex,parVortex,iG,iRc,ia,ib);
    if(err!=0){printf("problems\n"); return -5;}

    gsl_histogram_increment(hN,nCnect);

    err=histoIncVortex(nCnect,vCatalog,hG,hRc,ha,hb);
    if(err!=0){printf("problems\n"); return -5;}

    /* Preparing for printing */

    vortexQuickSort(parVortex,nVortex,&greaterAbsVorticity);
    vortexQuickSort(vCatalog,nCnect,&greaterAbsVorticity);

    err=fprintVortex(dadosVin,n,nVortex,parVortex);
    if(err!=0){printf("problems\n"); return -6;}

    err=fprintVortex(dadosVout,n,nCnect,vCatalog);
    if(err!=0){printf("problems\n"); return -6;}
  }

  fclose(dadosVin);
  fclose(dadosVout);

  sprintf(filename,"%s/histoOuG-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuRc-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOua-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOub-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuN-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);

  sprintf(filename,"%s/histoInG-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInRc-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoIna-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ia,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInb-%s.txt",folder,tag); dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ib,"%f","%f");
  fclose(dadosout);

  if(sField!=NULL)
    free(sField);
  if(gField!=NULL)
    free(gField);
  if(label!=NULL)
    free(label);
  if(vCatalog!=NULL)
    free(vCatalog);
  if(rCatalog!=NULL)
    free(rCatalog);
  if(majorVortex!=NULL)
    free(majorVortex);
  if(Glist!=NULL)
    free(Glist);
  if(Rclist!=NULL)
    free(Rclist);

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

  gsl_histogram_free(iG);
  gsl_histogram_free(iRc);
  gsl_histogram_free(ia);
  gsl_histogram_free(ib);
  /* histogram free - end*/

  return 0;
}

int fprintfRunParamSigned(FILE *dadosgen,long long int seed,float x0[],
                         float xf[],float dx[],float Gmin, float Gmax,
                         float rmin, float rmax, float xmin[], float xmax[], 
                         float v0y0)
{
  int i;

  if(dadosgen==NULL)
    return 1;

  fprintf(dadosgen,"seed: %d\n",seed);
  fprintf(dadosgen,"\ndomain: xi xf dx\n");
  fprintf(dadosgen,"%f %f %f\n",x0[0],xf[0],dx[0]);
  fprintf(dadosgen,"%f %f %f\n",x0[1],xf[1],dx[1]);
  fprintf(dadosgen,"\nvortex params: Uniform\n");
  fprintf(dadosgen,"G :+- %f %f\n",Gmin,Gmax);
  fprintf(dadosgen,"Rc: %f %f\n",rmin,rmax);
  fprintf(dadosgen,"a : %f %f\n",xmin[0],xmax[0]);
  fprintf(dadosgen,"b : %f %f\n",xmin[1],xmax[1]);
  fprintf(dadosgen,"\nshear v0/y0=%f\n",v0y0);
  
  return 0;
}

int histoIncVortex(int nVortex, float *parVortex,
                   gsl_histogram *iG, gsl_histogram *iRc,
                   gsl_histogram *ia, gsl_histogram *ib){
  int i;

  for(i=0;i<nVortex;i+=1){
    gsl_histogram_increment(iG,parVortex[4*i+0]);
    gsl_histogram_increment(iRc,parVortex[4*i+1]);
    gsl_histogram_increment(ia,parVortex[4*i+2]);
    gsl_histogram_increment(ib,parVortex[4*i+3]);
  }

  return 0;
}

int fprintVortex(FILE *dadosout, int run,int nVortex, float *vCatalog){
  int i;

  if(dadosout==NULL || run<0 || nVortex<=0 || vCatalog==NULL)
    return 1;

  for(i=0;i<nVortex;i+=1)
    fprintf(dadosout,"%d %d : %f %f %f %f\n",run,i,vCatalog[4*i+0]
                                                  ,vCatalog[4*i+1]
                                                  ,vCatalog[4*i+2]
                                                  ,vCatalog[4*i+3]);

  fprintf(dadosout,"\n");

  return 0;
}

int fprintsField(FILE *dadosout,float *x0,float *dx,
                 int Width, int Height, float *sField){
  int i,j;
  float x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || sField==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];

      fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
    }  

  return 0;
}

int fprintLabels(FILE *dadosout,float *x0,float *dx,
                 int Width, int Height, int *label){
  int i,j;
  float x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || label==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];

      fprintf(dadosout,"%f %f %d\n",x,y,label[i*Width+j]);
    }  

  return 0;
}