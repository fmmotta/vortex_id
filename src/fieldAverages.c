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
  double *gField=NULL,*uField=NULL,*uSubtr,*X,*Y,*wBkg,*background;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*ux=NULL,*uy=NULL;
  double *uAvg=NULL,*u2Avg=NULL,*eps=NULL,*epsSub=NULL,*epsPerp=NULL,*epsPerpSub=NULL;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1],bkgFile[400+1];
  FILE *dadosin,*dadosout;
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
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(uSubtr,2*Height*Width,double);
  fieldAlloc(background,2*Height*Width,double);
  fieldAlloc(  ux  ,2*Height*Width,double);
  fieldAlloc(  uy  ,2*Height*Width,double);
  fieldAlloc(wBkg,Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);
  fieldAlloc(uAvg,2*Height*Width,double);
  fieldAlloc(u2Avg,2*Height*Width,double);
  fieldAlloc(eps,Height*Width,double);
  fieldAlloc(epsSub,Height*Width,double);
  fieldAlloc(epsPerp,Height*Width,double);
  fieldAlloc(epsPerpSub,Height*Width,double);

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
    double x,y,Ux,Uy,sigmaUx,sigmaUy;
    double omega,strain,gamma,beta;
    FILE *dadosField;
    if(DEBUG_PRINT)
      printf("loading background file\n");
    dadosField=fopen(bkgFile,"r");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        fscanf(dadosField,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                                                          &x,&y,&Ux,&Uy,
                                                          &sigmaUx,&sigmaUy,
                                                          &omega,&gamma,
                                                          &beta,&strain);
        background[2*(i*Width+j)+0] = Ux;
        background[2*(i*Width+j)+1] = Uy;
        wBkg[i*Width+j] = omega;
      }

    fclose(dadosField);
  }

  dbgPrint(14,5);

  nSkip=0;

  //double *uAvg=NULL,*u2Avg=NULL,*eps=NULL,*epsSub=NULL,*epsPerp=NULL,*epsPerpSub=NULL;
  
  planeCount=0;
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uAvg[2*(i*Width+j)+0]=0.;
      uAvg[2*(i*Width+j)+1]=0.;
      u2Avg[2*(i*Width+j)+0]=0.;
      u2Avg[2*(i*Width+j)+1]=0.;
      eps[i*Width+j]=0.;
      epsSub[i*Width+j]=0.;
      epsPerp[i*Width+j]=0.;
      epsPerpSub[i*Width+j]=0.;
    }

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

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uAvg[2*(i*Width+j)+0] += uField[2*(i*Width+j)+0];
          uAvg[2*(i*Width+j)+1] += uField[2*(i*Width+j)+1];

          u2Avg[2*(i*Width+j)+0] += uSubtr[2*(i*Width+j)+0]*uSubtr[2*(i*Width+j)+0];
          u2Avg[2*(i*Width+j)+1] += uSubtr[2*(i*Width+j)+1]*uSubtr[2*(i*Width+j)+1];
        }

      dbgPrint(15,5);

      /****** Subtracted gradients and dissipation rate ******/
      
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

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          eps[i*Width+j] += ux[2*(i*Width+j)+0]*ux[2*(i*Width+j)+0] 
                         + uy[2*(i*Width+j)+1]*uy[2*(i*Width+j)+1]
                         + (ux[2*(i*Width+j)+0]+uy[2*(i*Width+j)+1])*(ux[2*(i*Width+j)+0]+uy[2*(i*Width+j)+1])
                         + (3./2.)*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1])*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1]);

          epsPerp[i*Width+j] += ux[2*(i*Width+j)+0]*ux[2*(i*Width+j)+0] 
                             + uy[2*(i*Width+j)+1]*uy[2*(i*Width+j)+1]
                             + (1./2.)*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1])*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1]);
        }      

      dbgPrint(15,7);
      /****** Subtracted gradients and dissipation rate ******/

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

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          epsSub[i*Width+j] += ux[2*(i*Width+j)+0]*ux[2*(i*Width+j)+0] 
                            + uy[2*(i*Width+j)+1]*uy[2*(i*Width+j)+1]
                            + (ux[2*(i*Width+j)+0]+uy[2*(i*Width+j)+1])*(ux[2*(i*Width+j)+0]+uy[2*(i*Width+j)+1])
                            + (3./2.)*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1])*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1]);

          epsPerpSub[i*Width+j] += ux[2*(i*Width+j)+0]*ux[2*(i*Width+j)+0] 
                                + uy[2*(i*Width+j)+1]*uy[2*(i*Width+j)+1]
                                + (1./2.)*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1])*(uy[2*(i*Width+j)+0]+ux[2*(i*Width+j)+1]);
        }      

      dbgPrint(15,9);  
    }
  } // End of Main loop

  printf("%d timesteps processed\n",n);

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uAvg[2*(i*Width+j)+0]  /= planeCount;
      uAvg[2*(i*Width+j)+1]  /= planeCount;
      u2Avg[2*(i*Width+j)+0] /= planeCount;
      u2Avg[2*(i*Width+j)+1] /= planeCount;
      eps[i*Width+j]         /= planeCount;
      epsSub[i*Width+j]      /= planeCount;
      epsPerp[i*Width+j]     /= planeCount;
      epsPerpSub[i*Width+j]  /= planeCount;
    }

  /***********************/

  printf("printing\n");

  double nu = 8.5764e-4;

  sprintf(filename,"%s/uAvg.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf %lf\n",X[j],Y[i],
                                           uAvg[2*(i*Width+j)+0],uAvg[2*(i*Width+j)+1]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/uBkg.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf\n",uAvg[2*(i*Width+j)+0],uAvg[2*(i*Width+j)+1]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/u2Avg.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf %lf\n",X[j],Y[i],u2Avg[2*(i*Width+j)+0],u2Avg[2*(i*Width+j)+1]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/eps.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf\n",X[j],Y[i],nu*eps[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsSub.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf\n",X[j],Y[i],nu*epsSub[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsPerp.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf\n",X[j],Y[i],nu*epsPerp[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsPerpSub.dat",folder);
  dadosout=fopen(filename,"w");
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      fprintf(dadosout,"%lf %lf %lf\n",X[j],Y[i],nu*epsPerpSub[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/uAvgLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf %lf\n",Y[i],uAvg[2*(i*Width+j)+0],uAvg[2*(i*Width+j)+1]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/u2AvgLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf %lf\n",Y[i],u2Avg[2*(i*Width+j)+0],u2Avg[2*(i*Width+j)+1]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf\n",Y[i],nu*eps[i*Width+j]);
  fclose(dadosout); dadosout=NULL;


  sprintf(filename,"%s/epsSubLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf\n",Y[i],nu*epsSub[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsPerpLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf\n",Y[i],nu*epsPerp[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/epsPerpSubLine.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1)
    fprintf(dadosout,"%lf %lf\n",Y[i],nu*epsPerpSub[i*Width+j]);
  fclose(dadosout); dadosout=NULL;

  sprintf(filename,"%s/kolmogorov.dat",folder);
  dadosout=fopen(filename,"w");
  j=0;
  for(i=0;i<Height;i+=1){
    double diss = nu*epsSub[i*Width+j];
    fprintf(dadosout,"%lf %lf %lf %lf ",Y[i],sqrt(sqrt(nu*nu*nu/diss)),sqrt(nu/diss),sqrt(sqrt(nu*diss)));
    diss = nu*eps[i*Width+j];
    fprintf(dadosout,"%lf %lf %lf \n",sqrt(sqrt(nu*nu*nu/diss)),sqrt(nu/diss),sqrt(sqrt(nu*diss)));
  }
  fclose(dadosout); dadosout=NULL;

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

  if(ux!=NULL) free(ux);
  if(uy!=NULL) free(uy);

  dbgPrint(23,0);

  if(wBkg!=NULL) free(wBkg);
  if(uAvg!=NULL) free(uAvg);
  if(u2Avg!=NULL) free(u2Avg);
  if(eps!=NULL) free(eps);
  if(epsSub!=NULL) free(epsSub);
  if(epsPerp!=NULL) free(epsPerp);
  if(epsPerpSub!=NULL) free(epsPerpSub);

  dbgPrint(24,0);

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
