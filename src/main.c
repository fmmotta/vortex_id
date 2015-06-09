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
  const int Width = 200, Height = 200,Pop=10,nVortex=20,nRuns=1000;
  const int numG=3,numRc=3;
  int seed=98755;
  int i,j,err,ngbr,found,nCnect,rCnect=0,*label,n,bin,nMax=500,pass=0;
  int nbList[8],eqList[Pop],**eqClass,it,nRecon=0;
  float Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  float xmin[2]={-9.,-9.},xmax[2]={9.,9.};
  float Glist[3]={1,5,10},Rclist[3]={0.5,1.0,1.5};
  float *parVortex=NULL,x0[2],dx[2],xf[2],*sField=NULL,*gField=NULL;
  float x,y,v0y0 = 0.00,*vCatalog=NULL,*rCatalog=NULL,*majorVortex=NULL;
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

  rCatalog = (float*)malloc(4*nMax*sizeof(float));
  if(rCatalog==NULL){
    printf("memory not allocked\n");
    return 4;
  }

  majorVortex = (float*)malloc(4*sizeof(float));
  if(majorVortex==NULL){
    printf("memory not allocked\n");
    return 5;
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

  dadosgen=fopen("data/multiRunGen.txt","w");
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
    err=genLOseenNaryList(numG,Glist,numRc,Rclist,xmin,xmax,seed,
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

    it=0;
    pass=0;
    rCnect=0;
    do{
      for(i=0;i<Height*Width;i+=1)
        label[i]=-1;

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

      vortexQuickSort(vCatalog,nCnect,&greaterCirculation);

      if(fabs(vCatalog[0])>threshold){
        pass=1; 
        if(rCnect>=nMax-1){
          printf("transpassing max size for rCatalog");
          return -10;
        }

        majorVortex[0]=-vCatalog[0]; rCatalog[4*rCnect+0] = vCatalog[0];
        majorVortex[1]= vCatalog[1]; rCatalog[4*rCnect+1] = vCatalog[1];
        majorVortex[2]= vCatalog[2]; rCatalog[4*rCnect+2] = vCatalog[2];
        majorVortex[3]= vCatalog[3]; rCatalog[4*rCnect+3] = vCatalog[3];

        rCnect+=1;
      }
      else
        break;

      err = addSingleOseen(1,majorVortex,x0,dx,Height,Width,&gField);
      if(err!=0){
        printf("alguma merda séria tá acontecendo\n");
        return err;
      }

      it+=1;

      for(i=0;i<4*nCnect;i+=1)
        vCatalog[i]=0.;
    }while( (pass!=0)&&(it<nVortex) ); //while(pass!=0);

    nRecon += rCnect;

    for(i=0;i<nVortex;i+=1){
      gsl_histogram_increment(iG,parVortex[4*i+0]);
      gsl_histogram_increment(iRc,parVortex[4*i+1]);
      gsl_histogram_increment(ia,parVortex[4*i+2]);
      gsl_histogram_increment(ib,parVortex[4*i+3]);
    }

    gsl_histogram_increment(hN,rCnect);
    for(i=0;i<rCnect;i+=1){
      gsl_histogram_increment(hG,rCatalog[4*i+0]);
      gsl_histogram_increment(hRc,rCatalog[4*i+1]);
      gsl_histogram_increment(ha,rCatalog[4*i+2]);
      gsl_histogram_increment(hb,rCatalog[4*i+3]);
    }
  }

  printf("nTotalVortex = %d nRecon=%d\n",nRuns*nVortex,nRecon);

  dadosout=fopen("data/histoOuG.txt","w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoOuRc.txt","w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoOua.txt","w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoOub.txt","w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoOuN.txt","w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);

  dadosout=fopen("data/histoInG.txt","w");
  gsl_histogram_fprintf(dadosout,iG,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoInRc.txt","w");
  gsl_histogram_fprintf(dadosout,iRc,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoIna.txt","w");
  gsl_histogram_fprintf(dadosout,ia,"%f","%f");
  fclose(dadosout);
  dadosout=fopen("data/histoInb.txt","w");
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
