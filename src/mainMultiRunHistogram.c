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
  const int Width = 200, Height = 200,Pop=10,nFixVortex=20,nRuns=100000;
  const int numG=3,numRc=3;
  int seed=98755,nVortex=20;
  int i,j,err,ngbr,found,nCnect,rCnect=0,*label,n,bin,nMax=500,pass=0;
  int nbList[8],eqList[Pop],**eqClass;
  float Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.5,threshold=0.5;
  float xmin[2]={-9.,-9.},xmax[2]={9.,9.};
  float Glist[3]={1,5,10},Rclist[3]={0.5,1.0,1.5};
  float *parVortex=NULL,x0[2],dx[2],xf[2],*sField=NULL,*gField=NULL;
  float x,y,v0y0 = 0.00,*vCatalog=NULL;
  FILE *dadosgen,*dadosout;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10;
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  gsl_histogram *iG,*iRc,*ia,*ib;
  
  threshold=0.5;

  x0[0]=-10.; xf[0]= 10.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=-10.; xf[1]= 10.; dx[1] = (xf[1]-x0[1])/Width;
  
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
  gsl_histogram_set_ranges_uniform(hG,-0.5,2.5*Gmax);
  hRc = gsl_histogram_alloc(hNRc); 
  gsl_histogram_set_ranges_uniform(hRc,-0.18,2.5*rmax);
  ha = gsl_histogram_alloc(hNa); 
  gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb = gsl_histogram_alloc(hNb); 
  gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN = gsl_histogram_alloc(hNN);
  gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);

  iG = gsl_histogram_alloc(hNG); 
  gsl_histogram_set_ranges_uniform(iG,0.,2.5*Gmax);
  iRc = gsl_histogram_alloc(hNRc); 
  gsl_histogram_set_ranges_uniform(iRc,0.,2.5*rmax);
  ia = gsl_histogram_alloc(hNa); 
  gsl_histogram_set_ranges_uniform(ia,xmin[0],xmax[0]);
  ib = gsl_histogram_alloc(hNb); 
  gsl_histogram_set_ranges_uniform(ib,xmin[1],xmax[1]);
  /* histogram preparation - end*/

  if(argc>1)
    seed = atoi(argv[1]);
  else
    seed = (int) time(NULL);

  seed = 98755; // Fix seed for comparison

  dadosgen=fopen("data/Uniform Comparison/Simple/multiRunGen.txt","w");
  fprintf(dadosgen,"seed: %d\n",seed);
  fprintf(dadosgen,"\ndomain: xi xf dx\n");
  fprintf(dadosgen,"%f %f %f\n",x0[0],xf[0],dx[0]);
  fprintf(dadosgen,"%f %f %f\n",x0[1],xf[1],dx[1]);
  fprintf(dadosgen,"\nvortex params: Binary\n");
  fprintf(dadosgen,"G : ");
  for(i=0;i<numG;i+=1)
    fprintf(dadosgen,"%f ",Glist[i]);
  fprintf(dadosgen,"\n");
  fprintf(dadosgen,"Rc: ");
  for(i=0;i<numG;i+=1)
    fprintf(dadosgen,"%f ",Rclist[i]);
  fprintf(dadosgen,"\n");
  fprintf(dadosgen,"a : %f %f\n",xmin[0],xmax[0]);
  fprintf(dadosgen,"b : %f %f\n",xmin[1],xmax[1]);
  fprintf(dadosgen,"\nshear v0/y0=%f\n",v0y0);
  fclose(dadosgen);

  for(n=0;n<nRuns;n+=1){
    if(n%1000 == 0)
      printf("%d runs have passed\n",n);
    //err=genLOseenNaryList(numG,Glist,numRc,Rclist,xmin,xmax,seed,
    //                      nVortex,&parVortex);
    nVortex = nFixVortex;
    err=genLOseenUniformList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
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

    for(i=0;i<nVortex;i+=1){
      gsl_histogram_increment(iG,parVortex[4*i+0]);
      gsl_histogram_increment(iRc,parVortex[4*i+1]);
      gsl_histogram_increment(ia,parVortex[4*i+2]);
      gsl_histogram_increment(ib,parVortex[4*i+3]);
    }

    gsl_histogram_increment(hN,nCnect);
    for(i=0;i<nCnect;i+=1){
      gsl_histogram_increment(hG,vCatalog[4*i+0]);
      gsl_histogram_increment(hRc,vCatalog[4*i+1]);
      gsl_histogram_increment(ha,vCatalog[4*i+2]);
      gsl_histogram_increment(hb,vCatalog[4*i+3]);
    }

  }

  dadosout=fopen("data/Uniform Comparison/Simple/histoOuG.txt","w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoOuRc.txt","w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoOua.txt","w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoOub.txt","w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoOuN.txt","w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);

  dadosout=fopen("data/Uniform Comparison/Simple/histoInG.txt","w");
  gsl_histogram_fprintf(dadosout,iG,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoInRc.txt","w");
  gsl_histogram_fprintf(dadosout,iRc,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoIna.txt","w");
  gsl_histogram_fprintf(dadosout,ia,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/Uniform Comparison/Simple/histoInb.txt","w");
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
