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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "ini.h"
#include "mt64.h"
#include "vortexGen.h"
#include "floodFill.h"
#include "stencilExtended.h"
#include "preprocessing.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "vortexExtractionExtend.h"
#include "inputManager.h"
#include "essayHandler.h"

#define DEBUG_MODE false
#define DEBUG_PRINT false
#define SUBTRACTION_MODE false

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

int readAxis(int Nx,int Ny,int Nz,int planeType,
             char *folder, double *X,double *Y);

int main(int argc,char **argv){
  //double cutoff; int rCnet=0;
  int Width = 100, Height = 100, Depth,nFixVortex=5;
  int runType=0,nSkip=0;
  int Nsnapshots=0,openFoamFile=0;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum=0,pln[8128],planeCount=0;
  int i,j,l,err,n,padWidth=2;//,k;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2];
  double t,t0,dt;
  double *datax,*datay;
  double *gField=NULL,*uField=NULL,*uSubtr,*X,*Y,*background;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*ux=NULL,*uy=NULL;
  double *uTilde,*uxTilde,*uyTilde,*yUAvg;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1],bkgFile[400+1];
  FILE *dadosin,*dadosout;
  gsl_fft_real_wavetable *realWTx,*realWTy;
  gsl_fft_real_workspace *realWSx,*realWSy;
  configVar cfg;

  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  err=initConfig(&cfg);
  
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0){
    printf("Can't load .ini file\n");
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
  planeIndex = cfg.pIndex;
  planeType  = cfg.pType;
  Nsnapshots = cfg.Nsnapshots;
  strcpy(foamFolder,cfg.FOAMfolder);
  
  if(planeIndex<0){
    if(cfg.planeNum>0){
      planeNum=cfg.planeNum;
      for(i=0;i<planeNum;i+=1)
        pln[i]=cfg.pln[i];
    }
    else{
      printf("Wrongly written configuration file, specify number of slices\n");
      return 1;
    }
  }else{
  	printf("This should not happen\n");
    return 4;
    planeNum = 1;
    pln[0] = planeIndex;
  }

       if(planeType==0){ Height = Ny; Width = Nx; Depth = Nz; }
  else if(planeType==1){ Height = Ny; Width = Nz; Depth = Nx; }
  else if(planeType==2){ Height = Nz; Width = Nx; Depth = Ny; }
  else{ printf("error, non-recognized plane type\n"); return -15; }

  if(planeIndex<0)
    printf("Switching to plane number list\n");

  if(planeIndex>=Depth)
    printf("Out of bounds plane\n");

  if(cfg.Nx == 0 || cfg.Ny == 0 || cfg.Nz == 0){
    printf("error, incompatible dimension sizes\n");
    return -16;
  }

  /* Loading Configuration -- I need something more concise */

  dbgPrint(2,0);

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
  
  runType = cfg.runType;
  
  dbgPrint(4,1);

  /**********************************/

  dbgPrint(7,0);

  if(runType==0)
    threshold = cfg.swThresh;
  else if (runType==2)
    threshold = cfg.sndSwThresh;
  else
    threshold = -1.;

  if(DEBUG_MODE)
    printf("threshold is %f\n",threshold);
  
  dbgPrint(8,0);

  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */

  dbgPrint(10,0);

  /* Memory Allocation */

  fieldAlloc(X,Nx+1,double);
  fieldAlloc(Y,Ny+1,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);

  fieldAlloc(gField ,4*Height*Width,double);
  fieldAlloc(uField ,2*Height*Width,double);
  fieldAlloc(uSubtr ,2*Height*Width,double);
  fieldAlloc(  ux   ,2*Height*Width,double);
  fieldAlloc(  uy   ,2*Height*Width,double);
  fieldAlloc(background,2*Height*Width,double);
  fieldAlloc(uBuff  ,2*(Height+2*padWidth)*(Width+2*padWidth),double);

  fieldAlloc(datax  ,Width ,double);
  fieldAlloc(datay  ,Height,double);
  fieldAlloc(yUAvg  ,2*Height,double);
  fieldAlloc(uTilde ,2*Height*Width,double);
  fieldAlloc(uxTilde,2*Height*Width,double);
  fieldAlloc(uyTilde,2*Height*Width,double);

  realWSx = gsl_fft_real_workspace_alloc (Width);
  realWTx = gsl_fft_real_wavetable_alloc (Width);
  realWSy = gsl_fft_real_workspace_alloc (Height);
  realWTy = gsl_fft_real_wavetable_alloc (Height);

  dbgPrint(13,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  dbgPrint(14,0);
  
  err=readAxis(Nx,Ny,Nz,planeType,folder,X,Y);

  dbgPrint(14,1);

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  dbgPrint(14,2);
  
  if(bkgFile[0]!='\0'){
    double Ux,Uy;
    FILE *dadosField;
    if(DEBUG_PRINT)
      printf("loading background file\n");
    dadosField=fopen(bkgFile,"r");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        fscanf(dadosField,"%lf %lf",&Ux,&Uy);
        background[2*(i*Width+j)+0] = Ux;
        background[2*(i*Width+j)+1] = Uy;
      }

    for(i=0;i<Height;i+=1){
      yUAvg[2*i+0] = 0.;
      yUAvg[2*i+1] = 0.;
      for(j=0;j<Width;j+=1){
        yUAvg[2*i+0] += background[2*(i*Width+j)+0];
        yUAvg[2*i+1] += background[2*(i*Width+j)+1];
      }
      yUAvg[2*i+0] /= Width;
      yUAvg[2*i+1] /= Width;
    }

    fclose(dadosField);
  }

  dbgPrint(14,5);

  nSkip=0;

  planeCount=0;

  for(n=0;n<Nsnapshots;n+=1){

    t=t0+((double)n)*dt;
    printf("%d timesteps processed\n",n);

    for(l=0;l<planeNum;l+=1){
      
      openFoamFile = 0;
      for(i=0;i<2*Height*Width;i+=1)
        uField[i]=0.;

      dbgPrint(15,0);

      if(DEBUG_PRINT)
        printf("plane =%d\n",pln[l]);

      if(planeType==0)      sprintf(filename,"%s/plane-z%d-%g.dat",folder,pln[l],t);
      else if(planeType==1) sprintf(filename,"%s/plane-x%d-%g.dat",folder,pln[l],t);
      else if(planeType==2) sprintf(filename,"%s/plane-y%d-%g.dat",folder,pln[l],t);
      
      dadosin=fopen(filename,"r");
      if(dadosin==NULL){
        openFoamFile=1;
        nSkip+=1;
        printf("Failed to open slice - %s\n",filename);
        break;
      }

      dbgPrint(15,31);

      planeCount+=1;
        
      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1)
          fscanf(dadosin,"%lf%lf",&(uField[2*(i*Width+j)+0]),&(uField[2*(i*Width+j)+1]));
      fclose(dadosin);
        
      if(openFoamFile!=0)
        continue;

      dbgPrint(15,3);

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uSubtr[2*(i*Width+j)+0]= uField[2*(i*Width+j)+0] 
                                 - background[2*(i*Width+j)+0];
          uSubtr[2*(i*Width+j)+1]= uField[2*(i*Width+j)+1]
                                 - background[2*(i*Width+j)+1];
        }

      dbgPrint(15,4);

      /*******************************************************/
      /***** UnSubtracted Field and Gradient Spectra *****/
      
      err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
      if(err!=0)
        return -2;
      
      err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
      if(err!=0)
        return -3;
      
      err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
      if(err!=0)
        return -4;
      
      dbgPrint(15,6);

      for(i=0;i<Height;i+=1){
        for(j=0;j<Width;j+=1)
          datax[j] = uField[2*(i*Width+j)+0] - yUAvg[2*i+0];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uTilde[2*(i*Width+j)+0] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];       

        for(j=0;j<Width;j+=1)
          datax[j] = uField[2*(i*Width+j)+1] - yUAvg[2*i+1];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uTilde[2*(i*Width+j)+1] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];

        for(j=0;j<Width;j+=1)
          datax[j] = ux[2*(i*Width+j)+0];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uxTilde[2*(i*Width+j)+0] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];       

        for(j=0;j<Width;j+=1)
          datax[j] = ux[2*(i*Width+j)+1];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uxTilde[2*(i*Width+j)+1] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];
        
        for(j=0;j<Width;j+=1)
          datax[j] = uy[2*(i*Width+j)+0];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uyTilde[2*(i*Width+j)+0] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];       

        for(j=0;j<Width;j+=1)
          datax[j] = uy[2*(i*Width+j)+1];
        gsl_fft_real_transform(datax,1,Width,realWTx, realWSx);
        for(j=0;j<(Width+1)/2;j+=1)
          uyTilde[2*(i*Width+j)+1] += datax[2*j]*datax[2*j]+datax[2*j+1]*datax[2*j+1];
      }


      /*******************************************************/
      /****** Subtracted Field and Gradient Spectra ******/
      
      dbgPrint(15,7);

      err = uFieldTouBuff(Height,Width,uSubtr,uBuff,padWidth);
      if(err!=0)
        return -2;
      
      err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
      if(err!=0)
        return -3;
      
      err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
      if(err!=0)
        return -4;

      dbgPrint(15,8);
    }
  } // End of Main loop

  for(i=0;i<Height;i+=1)
    for(j=0;j<(Width+1)/2;j+=1){
       uTilde[2*(i*Width+j)+0] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);
       uTilde[2*(i*Width+j)+1] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);

      uxTilde[2*(i*Width+j)+0] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);
      uxTilde[2*(i*Width+j)+1] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);

      uyTilde[2*(i*Width+j)+0] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);
      uyTilde[2*(i*Width+j)+1] /= Nsnapshots*planeNum*(X[Width-1]-X[0]);
    }

  printf("%d timesteps processed\n",n);

  /***********************/

  printf("printing\n");
  
  {
    int c;
    double norm=0.;
    for(i=0;i<Height;i+=1){
      norm = 0.;
      for(j=0;j<(Width+1)/2;j+=1)
        norm += uTilde[2*(i*Width+j)+0];

      for(j=0;j<(Width+1)/2;j+=1){
         uTilde[2*(i*Width+j)+0] /= norm;
         uTilde[2*(i*Width+j)+1] /= norm;
        uxTilde[2*(i*Width+j)+0] /= norm;
        uxTilde[2*(i*Width+j)+1] /= norm;
        uyTilde[2*(i*Width+j)+0] /= norm;
        uyTilde[2*(i*Width+j)+1] /= norm;
      }
    }
  }

  {
    double Re=1000, Yref=-1.0;
    char tildename[400+1];
    sprintf(tildename,"%s/spectra",folder);
    err=mkdir(tildename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(err!=0 && err!=-1){
      printf("error creating spectra directory - %d\n",err);
      return err;
    }
    
    FILE *uTildeFile; // = fopen("data/uTildeX.dat","w");
    
    for(i=0;i<Height;i+=1){
      sprintf(tildename,"%s/spectra/uTildeX-%.3f.dat",folder,Re*(Y[i]-Yref));
      uTildeFile = fopen(tildename,"w");
      for(j=0;j<(Width+1)/2;j+=1)
        fprintf(uTildeFile,"%lf %g %g %g %g %g %g\n",(2*M_PI*j)/(X[Width-1]-X[0])
                          , uTilde[2*(i*Width+j)+0], uTilde[2*(i*Width+j)+1]
                          ,uxTilde[2*(i*Width+j)+0],uxTilde[2*(i*Width+j)+1]
                          ,uyTilde[2*(i*Width+j)+0],uyTilde[2*(i*Width+j)+1]);
      fclose(uTildeFile);
    }

  }

  printf("printed\n");

  /***********************/
    
  dbgPrint(22,0);

  if(X!=NULL) free(X);
  if(Y!=NULL) free(Y);
  if(Xbuff!=NULL) free(Xbuff);
  if(Ybuff!=NULL) free(Ybuff);

  dbgPrint(22,1);

  if(uField!=NULL) free(uField);
  if(uSubtr!=NULL) free(uSubtr);
  if(uBuff!=NULL) free(uBuff);

  dbgPrint(22,2);

  if(background!=NULL) free(background);
  if(gField!=NULL)  free(gField);

  dbgPrint(22,3);

  if(datax!=NULL) free(datax);
  if(datay!=NULL) free(datay);
  if(ux!=NULL) free(ux);
  if(uy!=NULL) free(uy);
  if(uTilde!=NULL) free(uTilde);
  if(uxTilde!=NULL) free(uxTilde);
  if(uyTilde!=NULL) free(uyTilde);

  dbgPrint(23,0);

  return 0;
}

int readAxis(int Nx,int Ny,int Nz,int planeType,
             char *folder, double *X,double *Y)
{
  int i;
  char filename[400+1];
  FILE *nFile;

  if(planeType==0){
    sprintf(filename,"%s/Xaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Yaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }
  else if(planeType==1){
    sprintf(filename,"%s/Zaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Yaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }
  else if(planeType==2){
    sprintf(filename,"%s/Zaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Xaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }

  return 0;
}
