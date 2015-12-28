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
#include "stencilExtended.h"
#include "preprocessing.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "inputManager.h"
#include "essayHandler.h"

#define DEBUG_MODE false
#define DEBUG_PRINT false

#define dbgPrint(num,num2) if(DEBUG_PRINT) printf("check point - %d-%d\n",(num),(num2))

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));\
                                  if(ptr==NULL){                         \
                                    printf("memory not allocked\n");     \
                                    return 1;                            \
                                  }                                      \
                                  else{                                  \
                                    for(i=0;i<(size);i+=1)               \
                                      ptr[i]=(type) 0;                   \
                                  }                                      \

int main(int argc,char **argv){
  int Width = 100, Height = 100, Depth, nVortex=5,nFixVortex=5,nRuns=1000;
  int runType=0,genType=0,numG=3,numRc=3,*label=NULL,**eqClass=NULL;
  long long int seed=98755;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10,Nsnapshots,N;
  int Nu,Np,Nx,Ny,Nz,Nn,auHeight,auWidth,planeIndex,planeType;
  int i,j,k,err,nCnect,rCnect=0,n,it,nMax=1024,pass=0,padWidth=2;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2];
  double *parVortex=NULL,*Glist=NULL,*Rclist=NULL,cutoff=0.,t,t0,dt;
  double *sField=NULL,*gField=NULL,*g2Field=NULL,*uField=NULL;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*X=NULL,*Y=NULL,*Z=NULL,*ux=NULL;
  double *uy=NULL,*uxxy=NULL,*uxyy=NULL,*uxxx=NULL,*uyyy=NULL;
  double x,y,v0y0 = 0.00,*vCatalog=NULL,*rCatalog=NULL,*majorVortex=NULL;
  double hGmin=0.,hGmax=0.,hRcMin=0.,hRcMax=0.;
  char genFile[300+1],folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1];
  FILE *dadosgen,*dadosin,*dadosout,*dadosVin,*dadosVout,*uFile,*pFile,*nFile;
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  gsl_histogram *iG,*iRc,*ia,*ib;
  configVar cfg;
  openFoamIcoData *node=NULL;
  
  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  err=initConfig(&cfg);
  
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
  
  dbgPrint(2,0);

  if(DEBUG_MODE==true){
    err=printConfig(&cfg);
    if(err!=0)
      return err;
  }
  
  if(cfg.dim!=2){
    printf("Dimension is not 2 - Can't Follow\n");
    return 2;
  }
  
  /**** OpenFOAM parameters setting *****/
  
  Nx = cfg.Nx;
  Ny = cfg.Ny;
  Nz = cfg.Nz;
  t0 = cfg.t0;
  dt = cfg.dt;
  k = cfg.pIndex;
  planeType  = cfg.pType;
  Nsnapshots = cfg.Nsnapshots;
  strcpy(foamFolder,cfg.FOAMfolder);
  
  /*
  printf("%d %d %d %f %f %d %d %d\n",cfg.Nx,cfg.Ny,cfg.Nz,cfg.t0,cfg.dt
                                    ,cfg.pIndex,cfg.pType,cfg.Nsnapshots);*/

  if(cfg.pType==0){
    Height = Ny;
    Width  = Nx;
    Depth  = Nz;
  }
  else if(cfg.pType==1){
    Height = Ny;
    Width  = Nz;
    Depth  = Nx; 
  }
  else if(cfg.pType==2){
    Height = Nz;
    Width  = Nx;
    Depth  = Ny; 
  }
  else{
    printf("error, non-recognized plane type\n");
    return -15;
  }

  if(cfg.Nx == 0 || cfg.Ny == 0 || cfg.Nz == 0){
    printf("error, incompatible dimension sizes\n");
    return -16;
  }

  /* Loading Configuration -- I need something more concise */

  dbgPrint(2,0);

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }
  
  dbgPrint(3,0);

  //seed    = cfg.seed;
  //Width   = cfg.Width;
  //Height  = cfg.Height;
  //nRuns   = cfg.nRuns;
  runType = cfg.runType;
  //genType = cfg.genType;
  //nFixVortex = cfg.nVortex;

  //numG    = cfg.numG;
  //numRc   = cfg.numRc;
  
  dbgPrint(4,0);
  /*
  Glist   = (double*)malloc(numG*sizeof(double));
  if(Glist==NULL){printf("Can't allocate Glist\n"); return 3;}
  for(i=0;i<numG;i+=1)
    Glist[i] = cfg.Glist[i];

  Rclist   = (double*)malloc(numRc*sizeof(double));
  if(Rclist==NULL){printf("Can't allocate Glist\n"); return 3;}
  for(i=0;i<numRc;i+=1)
    Rclist[i]=cfg.Rclist[i];
  */

  dbgPrint(4,1);

  /**********************************/

  dbgPrint(7,0);
  /*
  Gmin    = cfg.Gmin;
  Gmax    = cfg.Gmax;
  rmin    = cfg.RcMin;
  rmax    = cfg.RcMax;*/
  v0y0    = cfg.v0y0;
  cutoff  = cfg.cutoff;

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
  
  dbgPrint(8,0);

  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */
  /* Memory Allocation */

  dbgPrint(10,0);

  fieldAlloc(label,Height*Width,int);
  fieldAlloc(sField ,Height*Width,double);
  fieldAlloc(gField ,4*Height*Width,double);
  fieldAlloc(g2Field,4*Height*Width,double);
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(  ux  ,2*Height*Width,double);
  fieldAlloc(  uy  ,2*Height*Width,double);
  fieldAlloc( uxxy ,2*Height*Width,double);
  fieldAlloc( uxyy ,2*Height*Width,double);
  fieldAlloc( uxxx ,2*Height*Width,double);
  fieldAlloc( uyyy ,2*Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);
  fieldAlloc(X,Nx+1,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Y,Ny+1,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);
  fieldAlloc(Z,Nz+1,double);

  dbgPrint(5,0);

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }

  vCatalog = (double*)malloc(4*nMax*sizeof(double));
  if(vCatalog==NULL){
    printf("memory not allocked\n");
    return 3;
  }
  for(i=0;i<4*nMax;i+=1)
    vCatalog[i]=-1.;
  
  rCatalog = (double*)malloc(4*nMax*sizeof(double));
  if(rCatalog==NULL){
    printf("memory not allocked\n");
    return 4;
  }
  for(i=0;i<4*nMax;i+=1)
    rCatalog[i]=-1.;
  
  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  dbgPrint(11,0);

  /* histogram preparation - begin */
  hG = gsl_histogram_alloc(hNG);   gsl_histogram_set_ranges_uniform(hG,hGmin,hGmax);
  hRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(hRc,hRcMin,hRcMax);
  ha = gsl_histogram_alloc(hNa);   gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb = gsl_histogram_alloc(hNb);   gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN = gsl_histogram_alloc(hNN);   gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);

  iG = gsl_histogram_alloc(hNG);   gsl_histogram_set_ranges_uniform(iG,hGmin,hGmax);
  iRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(iRc,hRcMin,hRcMax);
  ia = gsl_histogram_alloc(hNa);   gsl_histogram_set_ranges_uniform(ia,xmin[0],xmax[0]);
  ib = gsl_histogram_alloc(hNb);   gsl_histogram_set_ranges_uniform(ib,xmin[1],xmax[1]);
  /* histogram preparation - end*/

  dbgPrint(12,0);
  /*
  sprintf(genFile,"%s/genfile-%s.dat",folder,tag);
  dadosgen=fopen(genFile,"w");
  err=fprintfRunParamSigned(dadosgen,seed,x0,xf,dx,Gmin,Gmax,rmin,
                            rmax,xmin,xmax,v0y0);
  fclose(dadosgen);
  */

  //sprintf(filename,"%s/inputVortexes.txt",folder);
  //dadosVin = fopen(filename,"w");
  //sprintf(filename,"%s/outputVortexes.txt",folder);
  //dadosVout = fopen(filename,"w");

  dbgPrint(13,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  dbgPrint(14,0);

  sprintf(filename,"%s/constant/polyMesh/points",foamFolder);
  nFile = fopen(filename,"r");
  err=loadAxis(nFile,Nx,Ny,Nz,X,Y,Z);
  if(err!=0)
    return err;
  fclose(nFile);

  for(i=0;i<Height;i+=1)
    Y[i] = (Y[i]+Y[i+1])/2.;

  for(j=0;j<Width;j+=1)
    X[j] = (X[j]+X[j+1])/2.;

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  for(n=0;n<Nsnapshots;n+=1){
    
    t=t0+((double)n)*dt;

    if(n%1000 == 0){
      printf("%d runs have passed\n",n);
      //fflush(dadosVin);
      //fflush(dadosVout);
    }

    for(i=0;i<2*Height*Width;i+=1)
      uField[i]=0.;

    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;
    
    dbgPrint(15,0);
    
    sprintf(filename,"%s/%.4f/U",foamFolder,t);
    uFile = fopen(filename,"r");
    if(uFile==NULL)
      printf("problems opening uFile - %d\n",n);

    sprintf(filename,"%s/%.4f/p",foamFolder,t);
    pFile = fopen(filename,"r");
    if(pFile==NULL)
      printf("problems opening pFile - %d\n",n);
    
    dbgPrint(15,1);
    
    err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
    if(err!=0)
      printf("Problems with loadFields\n");
    
    //fclose(pFile); fclose(uFile);
    
    dbgPrint(15,2);

    dadosout=fopen("data/refU.dat","w");
    for(j=0;j<Height;j+=1)
      for(i=0;i<Width;i+=1){
        uField[2*(j*Width+i)+0] = node[id(i,j,k)].u;
        uField[2*(j*Width+i)+1] = node[id(i,j,k)].v;
        fprintf(dadosout,"%lf %lf\n",uField[2*(j*Width+i)+0]
                                  ,uField[2*(j*Width+i)+1]);
      }
    fclose(dadosout);
    
    dbgPrint(15,3);

    // change here
    //err=calcScalarField(runType,Height,Width,x0,dx,nVortex,parVortex,gField,v0y0,sField);
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sField);
    if(err!=0){
      printf("Error in calcScalarField\n");
      return err;
    }

    dbgPrint(15,4);

    err = floodFill(sField,Width,Height,eqClass,label);
    if(err!=0)
      printf("Problems in floodFill\n");

    err = renameLabels(Height,Width,label);
    if(err>0)
      nCnect=err;
    else
      printf("problems with renameLabels - %d\n",err);

    dbgPrint(15,4);

    if(n%1000==0){
      sprintf(filename,"%s/sField-%d.txt",folder,n);
      dadosout = fopen(filename,"w");
      fprintUsfield(dadosout,X,Y,Height,Width,sField);
      fclose(dadosout);

      sprintf(filename,"%s/labels-%d.txt",folder,n);
      dadosout = fopen(filename,"w");
      fprintUlabels(dadosout,X,Y,Height,Width,label);
      fclose(dadosout);
    }

    // change here
    //err=vortexReconstruction(runType,Height,Width,nCnect,x0,dx,sField,
    //                         gField,label,&vCatalog);
    err=vortexUReconstruction(runType,Height,Width,nCnect,X,Y,sField, 
                              gField,label,&vCatalog);
    if(err!=0){
      printf("problems in vortexReconstruction\n");
      return err;
    }
    
    sprintf(filename,"%s/vortices.dat",folder);
    dadosout=fopen(filename,"w");
    for(i=0;i<nCnect;i+=1)
      fprintf(dadosout,"%.12f %.12f %.8f %.8f\n",vCatalog[4*i+0],vCatalog[4*i+1]
                                            ,vCatalog[4*i+2],vCatalog[4*i+3]);
    fclose(dadosout);

    //vortexQuickSort(parVortex,nVortex,&greaterAbsCirculation);
    vortexQuickSort(vCatalog,nCnect,&greaterAbsCirculation);
    
    /* filtering by cutoff */
    /*
    rCnect=0;
    for(i=0;i<nCnect;i+=1){
      if(fabs(vCatalog[4*i+0])>cutoff){
        rCnect += 1;
        rCatalog[4*i+0]=vCatalog[4*i+0];
        rCatalog[4*i+1]=vCatalog[4*i+1];
        rCatalog[4*i+2]=vCatalog[4*i+2];
        rCatalog[4*i+3]=vCatalog[4*i+3];
      }
    }*/

    //err=histoIncVortex(nVortex,parVortex,iG,iRc,ia,ib);
    //if(err!=0){printf("problems\n"); return -5;}

    //gsl_histogram_increment(hN,nCnect);

    //err=histoIncVortex(rCnect,rCatalog,hG,hRc,ha,hb);
    //if(err!=0){printf("problems\n"); return -5;}
    dbgPrint(17,0);

    err=histoIncVortex(nCnect,vCatalog,hG,hRc,ha,hb);
    if(err!=0){printf("problems\n"); return -5;}
    
    dbgPrint(18,0);

    /* Preparing for printing */
    if(n%1000==0){
      dbgPrint(18,1);

      //sprintf(filename,"%s/vortexesIn-%.4f.txt",folder,t);
      //dadosVin = fopen(filename,"w");
      sprintf(filename,"%s/vortexesOu-%.4f.txt",folder,t);
      dadosVout = fopen(filename,"w");
      
      //err=fprintVortex(dadosVin,n,nVortex,parVortex);
      //if(err!=0){printf("problems fprintVortex\n"); return -6;}
  
      //err=fprintVortex(dadosVout,n,rCnect,rCatalog);
      //if(err!=0){printf("problems\n"); return -6;}
      
      dbgPrint(18,2);

      err=fprintVortex(dadosVout,n,nCnect,vCatalog);
      if(err!=0){printf("problems\n"); return -6;}
      
      dbgPrint(18,3);

      //fclose(dadosVin);
      fclose(dadosVout);
    }
  }
  
  dbgPrint(19,0);

  //fclose(dadosVin);
  //fclose(dadosVout);

  sprintf(filename,"%s/histoOuG-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuRc-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOua-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOub-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuN-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);
  
  dbgPrint(20,0);

  /*
  sprintf(filename,"%s/histoInG-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInRc-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoIna-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ia,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInb-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ib,"%f","%f");
  fclose(dadosout);
  */

  dbgPrint(21,0);

  sprintf(filename,"gnuplot_script.gnu");
  err=writeGnuplotScript(filename,folder,tag,nRuns,nVortex);
  if(err!=0){printf("Error printing gnuplot script\n");return err;}

  dbgPrint(22,0);

  if(sField!=NULL)
    free(sField);
  dbgPrint(22,1);
  if(gField!=NULL)
    free(gField);
  dbgPrint(22,2);
  if(label!=NULL)
    free(label);
  dbgPrint(22,3);
  if(vCatalog!=NULL)
    free(vCatalog);
  dbgPrint(22,4);
  //if(rCatalog!=NULL)
  //  free(rCatalog);
  if(ux!=NULL)
    free(ux);
  dbgPrint(22,5);
  if(uy!=NULL)
    free(uy);
  dbgPrint(22,6);
  if(uxxx!=NULL)
    free(uxxx);
  dbgPrint(22,7);
  if(uyyy!=NULL)
    free(uyyy);
  dbgPrint(22,8);
  if(uxxy!=NULL)
    free(uxxy);
  dbgPrint(22,9);
  if(uxyy!=NULL)
    free(uxyy);
  dbgPrint(22,10);
  //if(majorVortex!=NULL)
  //  free(majorVortex);
  //if(Glist!=NULL)
  //  free(Glist);
  //if(Rclist!=NULL)
  //  free(Rclist);
  
  dbgPrint(23,0);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);
  
  if(parVortex!=0)
    free(parVortex);

  dbgPrint(24,0);

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

/***********************/




/***********************/