#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mt64.h"
#include "floodFill.h"
#include "lambdaInit.h"
#include "vortexGen.h"
#include "vortexExtraction.h"

int main(int argc,char **argv){
  const int Width = 100, Height = 100, Pop=10,nVortex=5,nRuns=5000;
  int seed=98755;
  int i,j,err,ngbr,found,nCnect,*label,n;
  int nbList[8],eqList[Pop],**eqClass;
  float Gmin=1,Gmax=2,rmin=1,rmax=1;
  float xmin[2]={-9.,-9.},xmax[2]={9.,9.};
  float *parVortex=NULL,x0[2],dx[2],xf[2],*sField=NULL,*gField;
  float x,y,v0y0 = 0.00,*vCatalog=NULL;
  FILE *dadosout;
  
  dadosout=fopen("data/multiRunResults.txt","w");

  x0[0]=-10.; xf[0]= 10.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=-10.; xf[1]= 10.; dx[1] = (xf[1]-x0[1])/Width;

  if(argc>1)
    seed = atoi(argv[1]);
  else
    seed = (int) time(NULL);
 
  fprintf(dadosout,"seed: %d\n\n",seed);

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  gField = (float *)malloc(4*Height*Width*sizeof(float));
  if(gField==NULL){
    printf("memory not allocked\n");
    return 1;
  }
  
  label = (int*)malloc(Height*Width*sizeof(int));
  if(label==NULL){
    printf("memory not allocked\n");
    return 2;
  }
  
  //vCatalog = (float*)malloc(10*sizeof(float));

  for(n=0;n<nRuns;n+=1){
    err=genLOseenBinaryList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                             nVortex,&parVortex);
    if(err!=0)
      return err;

    for(i=0;i<4*Height*Width;i+=1)
      gField[i]=0.;

    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;

    err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0)
      printf("Problems in addSingleOseen\n");

    err=addConstXYShear(x0,dx,Height,Width,v0y0,&gField);
    if(err!=0)
      printf("Problems in addConstXYShear\n");

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

    fprintf(dadosout,"Input: G rc a b\n");
    for(i=0;i<nVortex;i+=1)
      fprintf(dadosout,"%lf %lf %lf %lf\n",parVortex[4*i+0],
                                           parVortex[4*i+1],
                                           parVortex[4*i+2],
                                           parVortex[4*i+3]);

    fprintf(dadosout,"Output: G rc a b\n");
    for(i=0;i<nCnect;i+=1)
      fprintf(dadosout,"%f %f %f %f\n",vCatalog[4*i+0],
                                       vCatalog[4*i+1],
                                       vCatalog[4*i+2],
                                       vCatalog[4*i+3]);
    fprintf(dadosout,"\n");
    free(vCatalog); vCatalog=NULL;
  }

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
  fclose(dadosout);

  return 0;
}
