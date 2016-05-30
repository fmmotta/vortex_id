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
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "vortexExtractionExtend.h"
#include "inputManager.h"
#include "essayHandler.h"

#define DEBUG_MODE false
#define DEBUG_PRINT false
#define PRINT_MODE false

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
  long long int seed=98755;
  int Width = 100, Height = 100,nVortex=5,nFixVortex=5,nRuns=1000;
  int runType=0,genType=0,numG=3,numRc=3,*label=NULL,**eqClass=NULL;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10, calcScalarMode=0,dataSize;
  int i,j,err,nCnect=0,rCnect=0,n,nMax=1024,padWidth=2,mCnect=0.;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2];
  double *parVortex=NULL,*Glist,*Rclist,cutoff=0.,*uAvgField,*u2AvgField;
  double *sField=NULL,*gField=NULL,*g2Field=NULL,*uField=NULL,*uSubtr=NULL;
  double *uBuff=NULL,*Xbuff,*Ybuff,*X,*Y,*mCatalog,*avgGradU,*background;
  double *ux,*uy,*uxxy,*uxyy,*uxxx,*uyyy,*sSubtr,*vortSndMomMatrix=NULL;
  double v0y0 = 0.00,*vCatalog=NULL,*rCatalog=NULL,*majorVortex=NULL,*wBkg=NULL;
  double hGmin=0.,hGmax=0.,hRcMin=0.,hRcMax=0.,sigmaUx,sigmaUy;
  char genFile[300+1],folder[100+1],tag[100+1],filename[400+1],bkgFile[400+1];
  FILE *dadosgen,*dadosout,*dadosVin,*dadosVout,*dadosField,*totalVin,*totalVout;
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  gsl_histogram *iG,*iRc,*ia,*ib;
  configVar cfg;

  dataSize = 12;
  
  /*********************************/

  dbgPrint(0,0);

  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  dbgPrint(1,0);

  /* Loading Configuration -- I need something more concise */

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
  
  dbgPrint(2,1);

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  if(cfg.bkgFile!=NULL)
    strcpy(bkgFile,cfg.bkgFile);
  else
    bkgFile[0]='\0';

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }
  
  dbgPrint(3,0);

  seed    = cfg.seed;
  Width   = cfg.Width;
  Height  = cfg.Height;
  nRuns   = cfg.nRuns;
  runType = cfg.runType;
  calcScalarMode = cfg.calcMode;
  genType = cfg.genType;
  nFixVortex = cfg.nVortex;

  numG    = cfg.numG;
  numRc   = cfg.numRc;
  
  dbgPrint(4,0);
  
  fieldAlloc(Glist,numG,double);
  fieldAlloc(Rclist,numRc,double);

  for(i=0;i<numG;i+=1)
    Glist[i] = cfg.Glist[i];

  for(i=0;i<numRc;i+=1)
    Rclist[i]=cfg.Rclist[i];
  
  dbgPrint(4,1);

  fieldAlloc(X,Width,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Y,Height,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);

  dbgPrint(5,0);

  /************************************/

  x0[0]   = cfg.x0[0]; x0[1] = cfg.x0[1]; 
  xmin[0] = x0[0]+1; xmin[1] = x0[1]+1;
  
  xf[0]   = cfg.xf[0]; xf[1] = cfg.xf[1]; 
  xmax[0] = xf[0]-1; xmax[1] = xf[1]-1;
  
  dx[0] = (xf[0]-x0[0])/Height;
  dx[1] = (xf[1]-x0[1])/Width;
  
  /***********************************/
  
  dbgPrint(6,0);

  for(j=0;j<Width;j+=1)
    X[j] = x0[0] + ((double)j)*dx[0];
  for(i=0;i<Height;i+=1)
    Y[i] = x0[1] + ((double)i)*dx[1];
  
  dbgPrint(6,1);

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  dbgPrint(6,2);

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  /**********************************/

  dbgPrint(7,0);

  Gmin    = cfg.Gmin;
  Gmax    = cfg.Gmax;
  rmin    = cfg.RcMin;
  rmax    = cfg.RcMax;
  v0y0    = cfg.v0y0;
  cutoff  = cfg.cutoff;

  if(runType==0)
    threshold = cfg.swThresh;
  else if (runType==2)
    threshold = cfg.sndSwThresh;
  else
    threshold = -1.;

  if(threshold<0)
    printf("negative threshold!\n");

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

  dbgPrint(9,0);

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  dbgPrint(10,0);
  
  fieldAlloc(vCatalog,4*nMax,double);
  fieldAlloc(rCatalog,4*nMax,double);
  fieldAlloc(mCatalog,dataSize*nMax,double);
  fieldAlloc(vortSndMomMatrix,4*nMax,double);
  fieldAlloc(avgGradU,4*nMax,double);
  fieldAlloc(sField ,Height*Width,double);
  fieldAlloc(sSubtr ,Height*Width,double);
  fieldAlloc(gField ,4*Height*Width,double);
  fieldAlloc(g2Field,4*Height*Width,double);
  fieldAlloc(label,Height*Width,int);
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(uSubtr,2*Height*Width,double);
  fieldAlloc(  ux  ,2*Height*Width,double);
  fieldAlloc(  uy  ,2*Height*Width,double);
  fieldAlloc( uxxy ,2*Height*Width,double);
  fieldAlloc( uxyy ,2*Height*Width,double);
  fieldAlloc( uxxx ,2*Height*Width,double);
  fieldAlloc( uyyy ,2*Height*Width,double);
  fieldAlloc(uAvgField,2*Height*Width,double);
  fieldAlloc(u2AvgField,2*Height*Width,double);
  fieldAlloc(background,2*Height*Width,double);
  fieldAlloc(wBkg,Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);

  dbgPrint(11,0);

  /* histogram preparation - begin */
  hG  = gsl_histogram_alloc(hNG);  gsl_histogram_set_ranges_uniform(hG,hGmin,hGmax);
  hRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(hRc,hRcMin,hRcMax);
  ha  = gsl_histogram_alloc(hNa);  gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb  = gsl_histogram_alloc(hNb);  gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN  = gsl_histogram_alloc(hNN);  gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);

  iG  = gsl_histogram_alloc(hNG);  gsl_histogram_set_ranges_uniform(iG,hGmin,hGmax);
  iRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(iRc,hRcMin,hRcMax);
  ia  = gsl_histogram_alloc(hNa);  gsl_histogram_set_ranges_uniform(ia,xmin[0],xmax[0]);
  ib  = gsl_histogram_alloc(hNb);  gsl_histogram_set_ranges_uniform(ib,xmin[1],xmax[1]);
  /* histogram preparation - end*/

  dbgPrint(12,0);

  sprintf(genFile,"%s/genfile-%s.dat",folder,tag);
  dadosgen=fopen(genFile,"w");
  err=fprintfRunParamSigned(dadosgen,seed,x0,xf,dx,Gmin,Gmax,rmin,
                            rmax,xmin,xmax,v0y0);
  fclose(dadosgen);

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

  if(DEBUG_PRINT){
    printf("v0y0=%lf\n",v0y0);
    printf("calcMode=%d\n",calcScalarMode);
  }

  for(i=0;i<2*Height*Width;i+=1)
    background[i]=0.;

  for(i=0;i<2*Height*Width;i+=1)
    uAvgField[i]=0.;

  for(i=0;i<2*Height*Width;i+=1)
    u2AvgField[i]=0.;
  
  if(bkgFile[0]!='\0'){
    double x,y,Ux,Uy;;
    double omega,strain,gamma,beta;
    if(DEBUG_PRINT)
      printf("loading background file\n");
    dadosField=fopen(bkgFile,"r");
    if(dadosField!=NULL){
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          fscanf(dadosField,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&Ux,&Uy,
                                                            &sigmaUx,&sigmaUy,
                                                            &omega,&gamma,
                                                            &beta,&strain);
          background[2*(i*Width+j)+0] = Ux;
          background[2*(i*Width+j)+1] = Uy;
          wBkg[i*Width+j] = omega;
        }
    }
    else{ // If no background file open, just zeroes the background field
      printf("No background file found\n");
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          background[2*(i*Width+j)+0] = Ux;
          background[2*(i*Width+j)+1] = Uy;
        }
    }
    if(dadosField!=NULL)
      fclose(dadosField);
  }

  sprintf(filename,"%s/inputVortexes.txt",folder);
  totalVin = fopen(filename,"w");
  if(totalVin==NULL){printf("Problems to open input vortex list\n");}

  sprintf(filename,"%s/outputVortexes.txt",folder);
  totalVout = fopen(filename,"w");
  if(totalVout==NULL){printf("Problems to open output vortex list\n");}

  for(n=0;n<nRuns;n+=1){

    if(n%1000 == 0){
      printf("%d runs have passed\n",n);
      //fflush(dadosVin);
      //fflush(dadosVout);
    }
    
    nVortex = nFixVortex;
    err=genVortices(genType,seed,xmin,xmax,nFixVortex,&parVortex,Gmin,Gmax,
                    rmin,rmax,numG,numRc,Glist,Rclist);
    if(err<0)
      return err;
    else if((err>0) && (err<nVortex))
      nVortex = err;
    
    if(calcScalarMode==0){
      if(DEBUG_PRINT)
        printf("Scalar Mode 0\n");
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uField[2*(i*Width+j)+0]= 0.;
          uField[2*(i*Width+j)+1]= 0.;
        }

      err=calcUScalarField(runType,Height,Width,padWidth,x0,dx,X,Y,Xbuff,
                           Ybuff,nVortex,parVortex,uField,uBuff,ux,uy,
                           uxxx,uyyy,uxxy,uxyy,gField,g2Field,
                           v0y0,sField);
    }
    else if(calcScalarMode==1){
      if(DEBUG_PRINT)
        printf("Scalar Mode 1\n");
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uField[2*(i*Width+j)+0]= 0.;
          uField[2*(i*Width+j)+1]= 0.;
        }

      err=calcScalarField(runType,Height,Width,x0,dx,nVortex,parVortex,
                          gField,v0y0,sField);
    }
    else if(calcScalarMode==2){
      if(DEBUG_PRINT)
        printf("Scalar Mode 2\n");
      
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uSubtr[2*(i*Width+j)+0]= 0. - background[2*(i*Width+j)+0];
          uSubtr[2*(i*Width+j)+1]= 0. - background[2*(i*Width+j)+1];
        }
      
      dbgPrint(15,0);
      err=calcUScalarField(runType,Height,Width,padWidth,x0,dx,X,Y,Xbuff,
                           Ybuff,nVortex,parVortex,uSubtr,uBuff,ux,uy,
                           uxxx,uyyy,uxxy,uxyy,gField,g2Field,
                           v0y0,sSubtr);
      if(err!=0)
        break;

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uField[2*(i*Width+j)+0]= 0.;
          uField[2*(i*Width+j)+1]= 0.;
        }

      dbgPrint(15,1);
      err=calcUScalarField(runType,Height,Width,padWidth,x0,dx,X,Y,Xbuff,
                           Ybuff,nVortex,parVortex,uField,uBuff,ux,uy,
                           uxxx,uyyy,uxxy,uxyy,gField,g2Field,
                           v0y0,sField);
      if(err!=0)
        break;

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          if( sField[i*Width+j] > sSubtr[i*Width+j])
            sField[i*Width+j] = sSubtr[i*Width+j];

          uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
          uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
        }
      dbgPrint(15,2);
    }
    else if(calcScalarMode==3){
      if(DEBUG_PRINT)
        printf("Scalar Mode 3\n");

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uField[2*(i*Width+j)+0]= 0.;
          uField[2*(i*Width+j)+1]= 0.;
        }

      err=calcUScalarField(runType,Height,Width,padWidth,x0,dx,X,Y,Xbuff,
                           Ybuff,nVortex,parVortex,uField,uBuff,ux,uy,
                           uxxx,uyyy,uxxy,uxyy,gField,g2Field,
                           v0y0,sField);
      if(err!=0)
        break;
      
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uSubtr[2*(i*Width+j)+0]= 0. - background[2*(i*Width+j)+0];
          uSubtr[2*(i*Width+j)+1]= 0. - background[2*(i*Width+j)+1];
        }

      err=calcUScalarField(runType,Height,Width,padWidth,x0,dx,X,Y,Xbuff,
                           Ybuff,nVortex,parVortex,uSubtr,uBuff,ux,uy,
                           uxxx,uyyy,uxxy,uxyy,gField,g2Field,
                           v0y0,sSubtr);
      if(err!=0)
        break;

      // gField have already been updated to the subtracted field. 
      // __ Do not remove this message __

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          if( sField[i*Width+j] > sSubtr[i*Width+j])
            sField[i*Width+j] = sSubtr[i*Width+j];

          uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
          uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
        }
      
    }
    else{
      printf("Not identified operation mode - %d\n",calcScalarMode);
      return -18;
    }

    if(err!=0){
      printf("Error in calcScalarField\n");
      return err;
    }

    dbgPrint(15,3);
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        uAvgField[2*(i*Width+j)+0] += uField[2*(i*Width+j)+0];
        uAvgField[2*(i*Width+j)+1] += uField[2*(i*Width+j)+1];

        u2AvgField[2*(i*Width+j)+0] += uField[2*(i*Width+j)+0]*uField[2*(i*Width+j)+0];
        u2AvgField[2*(i*Width+j)+1] += uField[2*(i*Width+j)+1]*uField[2*(i*Width+j)+1];
      }

    dbgPrint(15,4);
    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;

    dbgPrint(15,5);
    err = floodFill(sField,Width,Height,eqClass,label);
    if(err!=0)
      printf("Problems in floodFill\n");
    
    err = renameLabels(Height,Width,label);
    if(err>0)
      nCnect=err;
    else
      printf("problems with renameLabels - %d\n",err);
    
    dbgPrint(15,6);

    if((n%1000==0) && PRINT_MODE){
      sprintf(filename,"%s/sField-%s-%d.txt",folder,tag,n);
      dadosField = fopen(filename,"w");
      fprintsField(dadosField,x0,dx,Height,Width,sField);
      fclose(dadosField);

      sprintf(filename,"%s/labels-%s-%d.txt",folder,tag,n);
      dadosField = fopen(filename,"w");
      fprintLabels(dadosField,x0,dx,Width,Height,label);
      fclose(dadosField);

    }

    dbgPrint(15,7);
    // WARNING : Change Here
    //err=vortexReconstruction(runType,Height,Width,nCnect,x0,dx,sField,
    //                         gField,label,&vCatalog);
    
    
    //err=vortexUReconstruction(runType,Height,Width,nCnect,X,Y,sField, 
    //                          gField,label,&vCatalog);
    //if(err!=0){
    //  printf("problems in vortexReconstruction\n");
    //  return err;
    //}
    
    err=extract012Momentsw2(Height,Width,nCnect,X,Y,sField,gField,label,
                            vCatalog,vortSndMomMatrix,avgGradU);
    if(err!=0){
      printf("problems in extract012Momentsw2\n");
      return err;
    }
    
    dbgPrint(15,8);
    if(calcScalarMode==2){
      double bkgG[nCnect];
      
      for(i=0;i<nCnect;i+=1)
        bkgG[i] = 0.;

      err = extractAvgBkgVort(Height,Width,X,Y,nCnect,label,wBkg,bkgG);
      for(i=0;i<nCnect;i+=1)
        vCatalog[4*i+0] = vCatalog[4*i+0]-bkgG[i];
    }
    dbgPrint(15,9);
    
    // Here for vortex position discart

    for(i=0;i<nCnect;i+=1){
      //vCatalog[4*i+0] += M_PI*vCatalog[4*i+1]*vCatalog[4*i+1]*v0y0;
      if(runType==0){
        vCatalog[4*i+0]= 1.397948086*vCatalog[4*i+0];
        vCatalog[4*i+1]=(1./1.12091)*vCatalog[4*i+1];
      }
      else if(runType==1){
        vCatalog[4*i+0]= 2.541494083*vCatalog[4*i+0];
        vCatalog[4*i+1]=  (sqrt(2.))*vCatalog[4*i+1]; 
      }
    }

    for(i=0;i<nCnect;i+=1){
      mCatalog[dataSize*i+0]  = vCatalog[4*i+0];
      mCatalog[dataSize*i+1]  = vCatalog[4*i+1];
      mCatalog[dataSize*i+2]  = vCatalog[4*i+2];
      mCatalog[dataSize*i+3]  = vCatalog[4*i+3];
      mCatalog[dataSize*i+4]  = vortSndMomMatrix[4*i+0];
      mCatalog[dataSize*i+5]  = vortSndMomMatrix[4*i+1];
      mCatalog[dataSize*i+6]  = vortSndMomMatrix[4*i+2];
      mCatalog[dataSize*i+7]  = vortSndMomMatrix[4*i+3];
      mCatalog[dataSize*i+8]  = avgGradU[4*i+0];
      mCatalog[dataSize*i+9]  = avgGradU[4*i+1];
      mCatalog[dataSize*i+10] = avgGradU[4*i+2];
      mCatalog[dataSize*i+11] = avgGradU[4*i+3];
    }

    vortexQuickSort(parVortex,nVortex,&greaterAbsCirculation);
    vortexQuickSort(vCatalog,nCnect,&greaterAbsCirculation);
    vortexAdaptiveQuickSort(mCatalog,nCnect,dataSize,&greaterAbsCirculation);

    /* filtering by cutoff */

    rCnect=0;
    for(i=0;i<nCnect;i+=1){
      if(fabs(vCatalog[4*i+0])>cutoff){
        rCnect += 1;
        rCatalog[4*i+0]=vCatalog[4*i+0];
        rCatalog[4*i+1]=vCatalog[4*i+1];
        rCatalog[4*i+2]=vCatalog[4*i+2];
        rCatalog[4*i+3]=vCatalog[4*i+3];
      }
    }

    mCnect=0;
    for(i=0;i<nCnect;i+=1)
      if(fabs(mCatalog[dataSize*i+0])>cutoff)
        mCnect += 1;

    /* printing to histogram */

    err=histoIncVortex(nVortex,parVortex,iG,iRc,ia,ib);
    if(err!=0){printf("problems\n"); return -5;}

    gsl_histogram_increment(hN,rCnect);

    err=histoIncVortex(rCnect,rCatalog,hG,hRc,ha,hb);
    if(err!=0){printf("problems\n"); return -5;}

    /* Preparing for printing */
    if((n%1000==0) && PRINT_MODE){
      sprintf(filename,"%s/vortexesIn-%d.txt",folder,n);
      dadosVin = fopen(filename,"w");
      sprintf(filename,"%s/vortexesOu-%d.txt",folder,n);
      dadosVout = fopen(filename,"w");
      
      err=fprintVortex(dadosVin,n,nVortex,parVortex);
      if(err!=0){printf("problems\n"); return -6;}
  
      err=fprintVortex(dadosVout,n,rCnect,rCatalog);
      if(err!=0){printf("problems\n"); return -6;}

      fclose(dadosVin);
      fclose(dadosVout);

      sprintf(filename,"%s/vortexMoments-%d.txt",folder,n);
      dadosVout = fopen(filename,"w");
      err=fprintSafeVortexMoments(dadosVout,n,dataSize,mCnect,mCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems vortexSafeMoments\n"); return -6;}
      fclose(dadosVout);
    }

    {
      fprintVortex(totalVin,n,nVortex,parVortex);
      fprintVortex(totalVout,n,rCnect,rCatalog);

      if(n%1000==0){
        fflush(totalVin);
        fflush(totalVout);
      }
    }
  }

  if(totalVin!=NULL) fclose(totalVin);
  if(totalVout!=NULL) fclose(totalVout);

  printf("%d runs have passed\n",n);
  
  {
    double omega,strain,gamma,beta;

    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        uAvgField[2*(i*Width+j)+0] /= nRuns;
        uAvgField[2*(i*Width+j)+1] /= nRuns;

        u2AvgField[2*(i*Width+j)+0] /= nRuns;
        u2AvgField[2*(i*Width+j)+1] /= nRuns;
      }

    err = uFieldTouBuff(Height,Width,uAvgField,uBuff,padWidth);
    if(err!=0)
      return -2;
  
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -3;

    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -4;

    sprintf(filename,"%s/uAvgField-%s.txt",folder,tag); 
    dadosout = fopen(filename,"w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        fprintf(dadosout,"%f %f %f %f ",X[j],Y[i]
                                       ,uAvgField[2*(i*Width+j)+0]
                                       ,uAvgField[2*(i*Width+j)+1]);
        
        sigmaUx = u2AvgField[2*(i*Width+j)+0] - 
                  uAvgField[2*(i*Width+j)+0]*uAvgField[2*(i*Width+j)+0];
        
        sigmaUy = u2AvgField[2*(i*Width+j)+1] - 
                  uAvgField[2*(i*Width+j)+1]*uAvgField[2*(i*Width+j)+1];
        fprintf(dadosout,"%f %f ",sigmaUx,sigmaUy);
        
        omega = ux[2*(i*Width+j)+1]-uy[2*(i*Width+j)+0];
        gamma = ux[2*(i*Width+j)+1]+uy[2*(i*Width+j)+0];
        beta  = ux[2*(i*Width+j)+0];
        strain = sqrt(gamma*gamma+beta*beta);

        fprintf(dadosout,"%f %f %f %f",omega,gamma,beta,strain);

        fprintf(dadosout,"\n");
      }
    fclose(dadosout);

    sprintf(filename,"%s/vertical-x0-%s.txt",folder,tag);
    dadosout = fopen(filename,"w");
    j=Width/2;
    for(i=0;i<Height;i+=1){
      fprintf(dadosout,"%f %f %f\n",Y[i],uAvgField[2*(i*Width+j)+0]
                                        ,uAvgField[2*(i*Width+j)+1]);

      sigmaUx = u2AvgField[2*(i*Width+j)+0] - 
                uAvgField[2*(i*Width+j)+0]*uAvgField[2*(i*Width+j)+0];

      sigmaUy = u2AvgField[2*(i*Width+j)+1] - 
                uAvgField[2*(i*Width+j)+1]*uAvgField[2*(i*Width+j)+1];
      fprintf(dadosout,"%f %f \n",sigmaUx,sigmaUy);
    }
    fclose(dadosout);
  }

  dbgPrint(16,0);

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

  sprintf(filename,"gnuplot_script.gnu");
  err=writeGnuplotScript(filename,folder,tag,nRuns,nVortex);
  if(err!=0){printf("Error printing gnuplot script\n");return err;}
  
  dbgPrint(17,0);

  if(X!=NULL) free(X);
  if(Y!=NULL) free(Y);
  if(Xbuff!=NULL) free(Xbuff);
  if(Ybuff!=NULL) free(Ybuff);
  if(ux!=NULL) free(ux);
  if(uy!=NULL) free(uy);
  if(uxxx!=NULL) free(uxxx);
  if(uxyy!=NULL) free(uxyy);
  if(uxxy!=NULL) free(uxxy);
  if(uyyy!=NULL) free(uyyy);
  if(uField!=NULL) free(uField);
  if(uSubtr!=NULL) free(uSubtr);
  if(uBuff!=NULL) free(uBuff);
  if(uAvgField!=NULL) free(uAvgField);
  if(u2AvgField!=NULL) free(u2AvgField);
  if(background!=NULL) free(background);
  if(sField!=NULL) free(sField);
  if(sSubtr!=NULL) free(sSubtr);
  if(wBkg!=NULL) free(wBkg);
  if(gField!=NULL) free(gField);
  if(g2Field!=NULL) free(g2Field);
  if(label!=NULL)  free(label);
  if(vCatalog!=NULL) free(vCatalog);
  if(rCatalog!=NULL) free(rCatalog);
  if(mCatalog!=NULL) free(mCatalog);
  if(vortSndMomMatrix!=NULL) free(vortSndMomMatrix);
  if(avgGradU!=NULL) free(avgGradU);
  if(majorVortex!=NULL) free(majorVortex);
  if(Glist!=NULL) free(Glist);
  if(Rclist!=NULL) free(Rclist);

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


    /*
    {
      rCnect=0;
      for(i=0;i<nCnect;i+=1){
        if(fabs(vCatalog[4*i+2])<=9 && fabs(vCatalog[4*i+3])<=9){
          rCatalog[4*rCnect+0]=vCatalog[4*i+0];
          rCatalog[4*rCnect+1]=vCatalog[4*i+1];
          rCatalog[4*rCnect+2]=vCatalog[4*i+2];
          rCatalog[4*rCnect+3]=vCatalog[4*i+3];
          rCnect+=1;
        }
      }

      for(i=0;i<rCnect;i+=1){
        vCatalog[4*i+0]=rCatalog[4*i+0];
        vCatalog[4*i+1]=rCatalog[4*i+1];
        vCatalog[4*i+2]=rCatalog[4*i+2];
        vCatalog[4*i+3]=rCatalog[4*i+3];      
      }
      for(i=rCnect;i<nCnect;i+=1){
        vCatalog[4*i+0]=0.;
        vCatalog[4*i+1]=0.;
        vCatalog[4*i+2]=0.;
        vCatalog[4*i+3]=0.;      
      }
      nCnect=rCnect;
    }
    */