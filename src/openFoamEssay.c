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

#define DEBUG_MODE true
#define DEBUG_PRINT true

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
  int runType=0,*label=NULL,**eqClass=NULL;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10,Nsnapshots;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum,pln[8128];
  int i,j,k,l,err,nCnect=0,rCnect=0,n,nMax=1024,padWidth=2;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2];
  double *parVortex=NULL,cutoff=0.,t,t0,dt;
  double *sField=NULL,*gField=NULL,*g2Field=NULL,*uField=NULL,*X,*Y;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*Xload=NULL,*Yload=NULL;
  double *Zload=NULL,*ux=NULL;
  double *uy=NULL,*uxxy=NULL,*uxyy=NULL,*uxxx=NULL,*uyyy=NULL;
  double v0y0 = 0.00,*vCatalog=NULL,*rCatalog=NULL;
  double hGmin=0.,hGmax=0.,hRcMin=0.,hRcMax=0.;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1];
  FILE *dadosout,*dadosVout,*uFile,*pFile,*nFile,*vortexFile;
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
  }

  if(planeType==0){
    Height = Ny;
    Width  = Nx;
    Depth  = Nz;
  }
  else if(planeType==1){
    Height = Ny;
    Width  = Nz;
    Depth  = Nx; 
  }
  else if(planeType==2){
    Height = Nz;
    Width  = Nx;
    Depth  = Ny; 
  }
  else{
    printf("error, non-recognized plane type\n");
    return -15;
  }

  if(planeIndex<0 || planeIndex>=Depth)
    printf("Out of bounds plane\n");

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
  
  runType = cfg.runType;
  
  dbgPrint(4,1);

  /**********************************/

  dbgPrint(7,0);
  
  v0y0    = cfg.v0y0;
  cutoff  = cfg.cutoff;

  if(runType==0)
    threshold = cfg.swThresh;
  else if (runType==2)
    threshold = cfg.sndSwThresh;
  else
    threshold = -1.;

  if(DEBUG_MODE)
    printf("threshold is %f\n",threshold);

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
  fieldAlloc(Y,Ny+1,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);
  fieldAlloc(Xload,Nx+1,double);
  fieldAlloc(Yload,Ny+1,double);
  fieldAlloc(Zload,Nz+1,double);

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

  dbgPrint(13,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  dbgPrint(14,0);
  
  if(planeIndex>=0){
    sprintf(filename,"%s/constant/polyMesh/points",foamFolder);
    nFile = fopen(filename,"r");
    err=loadAxis(nFile,Nx,Ny,Nz,Xload,Yload,Zload);
    if(err!=0)
      return err;
    fclose(nFile);

    if(planeType==0){
      for(i=0;i<Height;i+=1)
        Y[i] = (Yload[i]+Yload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Xload[j]+Xload[j+1])/2.;
    }
    else if(planeType==1){
      for(i=0;i<Height;i+=1)
        Y[i] = (Yload[i]+Yload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Zload[j]+Zload[j+1])/2.; 
    }
    else if(planeType==2){
      for(i=0;i<Height;i+=1)
        Y[i] = (Zload[i]+Zload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Xload[j]+Xload[j+1])/2.; 
    }
    else 
      printf("non-identified plane type\n");
  }
  else{
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
      sprintf(filename,"%s/Yaxis.dat",folder);
      nFile=fopen(filename,"r");
      for(i=0;i<Nx;i+=1)
        fscanf(nFile,"%lf",&(X[i]));
      fclose(nFile); nFile=NULL;

      sprintf(filename,"%s/Zaxis.dat",folder);
      nFile=fopen(filename,"r");
      for(i=0;i<Ny;i+=1)
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
  }

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");
  
  sprintf(filename,"%s/vortices.dat",folder);
  vortexFile = fopen(filename,"w");
  
  for(n=0;n<Nsnapshots;n+=1){
    
    t=t0+((double)n)*dt;

    //if(n%10 == 0){
      printf("%d runs have passed\n",n);
      fflush(vortexFile);
    //}

    for(i=0;i<2*Height*Width;i+=1)
      uField[i]=0.;

    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;
    
    dbgPrint(15,0);
    if(planeIndex>=0){
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
    
      fclose(pFile); fclose(uFile);
    
      dbgPrint(15,2);
    
      if(DEBUG_PRINT)
        printf("folder = %s\n",folder);
      sprintf(filename,"%s/refU-%.4f.dat",folder,t);
      dadosout=NULL;
      dadosout=fopen(filename,"w");
      if(dadosout==NULL){
        printf("Could not open Folder\n");
        perror("Error opening reference file\n");
        exit(EXIT_FAILURE);
      }

      if(planeType==0){
        if(DEBUG_PRINT)
          printf("XY plane\n");
        k=planeIndex;
        for(j=0;j<Height;j+=1)
          for(i=0;i<Width;i+=1){
            uField[2*(j*Width+i)+0] = node[id(i,j,k)].u;
            uField[2*(j*Width+i)+1] = node[id(i,j,k)].v;
            fprintf(dadosout,"%lf %lf\n",uField[2*(j*Width+i)+0]
                                        ,uField[2*(j*Width+i)+1]);
          }
      }
      else if(planeType==1){
        if(DEBUG_PRINT)
          printf("YZ plane\n");
        i=planeIndex;
        for(j=0;j<Height;j+=1)
          for(k=0;k<Width;k+=1){
            uField[2*(j*Width+k)+0] = node[id(i,j,k)].w;
            uField[2*(j*Width+k)+1] = node[id(i,j,k)].v;
            fprintf(dadosout,"%lf %lf\n",uField[2*(j*Width+k)+0]
                                        ,uField[2*(j*Width+k)+1]);
          }
      }
      else if(planeType==2){
        if(DEBUG_PRINT)
          printf("ZX plane\n");
        j=planeIndex;
        for(k=0;k<Height;k+=1)
          for(i=0;i<Width;i+=1){
            uField[2*(k*Width+i)+0] = node[id(i,j,k)].w;
            uField[2*(k*Width+i)+1] = node[id(i,j,k)].u;
            fprintf(dadosout,"%lf %lf\n",uField[2*(k*Width+i)+0]
                                        ,uField[2*(k*Width+i)+1]);
          }
      }
    
      if(dadosout!=NULL)
        fclose(dadosout);
    
    }
    else{
      // WARNING - This section is not ready to use!!!
      if(planeType==0){
        if(DEBUG_PRINT)
          printf("XY plane\n");
        for(l=0;l<planeNum;l+=1){
          k=pln[l];
          sprintf(filename,"%s/plane-z%3d-%.4f.dat",folder,k,t);

        }
      }
    }

    dbgPrint(15,3);

    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sField);
    if(err!=0){
      printf("Error in calcScalarField - %d\n",err);
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

    if(n%10==0){
      sprintf(filename,"%s/sField-%d.txt",folder,n);
      dadosout = fopen(filename,"w");
      fprintUsfield(dadosout,X,Y,Height,Width,sField);
      fclose(dadosout);

      sprintf(filename,"%s/labels-%d.txt",folder,n);
      dadosout = fopen(filename,"w");
      fprintUlabels(dadosout,X,Y,Height,Width,label);
      fclose(dadosout);

      sprintf(filename,"%s/presence-%d.txt",folder,n);
      dadosout = fopen(filename,"w");
      fprintUpresence(dadosout,X,Y,Height,Width,label);
      fclose(dadosout);
    }

    err=vortexUReconstruction(runType,Height,Width,nCnect,X,Y,sField, 
                              gField,label,&vCatalog);
    if(err!=0){
      printf("problems in vortexReconstruction\n");
      return err;
    }
    /*
    sprintf(filename,"%s/vortices-%d.dat",folder,n);
    dadosout=fopen(filename,"w");
    for(i=0;i<nCnect;i+=1)
      fprintf(dadosout,"%.12f %.12f %.8f %.8f\n",vCatalog[4*i+0],vCatalog[4*i+1]
                                            ,vCatalog[4*i+2],vCatalog[4*i+3]);
    fclose(dadosout);*/

    vortexQuickSort(vCatalog,nCnect,&greaterAbsCirculation);
    
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

    dbgPrint(17,0);

    err=histoIncVortex(nCnect,vCatalog,hG,hRc,ha,hb);
    if(err!=0){printf("problems\n"); return -5;}
    
    dbgPrint(18,0);

    /* Preparing for printing */
    if(n%10==0){

      sprintf(filename,"%s/vortices-%.4f.txt",folder,t);
      dadosVout = fopen(filename,"w");   
      err=fprintVortex(dadosVout,n,nCnect,vCatalog);
      if(err!=0){printf("problems\n"); return -6;}
      fclose(dadosVout);

      
      sprintf(filename,"%s/vorticesSafe-%.4f.txt",folder,t);
      dadosVout = fopen(filename,"w");
      err=fprintSafeVortex(dadosVout,n,nCnect,vCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems\n"); return -6;}
      fclose(dadosVout);

    }

    err=fprintSafeVortex(vortexFile,n,nCnect,vCatalog,Height,Width,X,Y);
    if(err!=0){printf("problems in printing vortexfile\n"); return -6;}
  } // End of Main loop

  fclose(vortexFile);
  
  dbgPrint(19,0);

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

  sprintf(filename,"gnuplot_script.gnu");
  err=writeGnuplotScript(filename,folder,tag,nRuns,nVortex);
  if(err!=0){printf("Error printing gnuplot script\n");return err;}

  dbgPrint(22,0);

  if(X!=NULL) free(X);
  if(Y!=NULL) free(Y);
  if(Xbuff!=NULL) free(Xbuff);
  if(Ybuff!=NULL) free(Ybuff);
  if(sField!=NULL) free(sField); 
  if(gField!=NULL)  free(gField);
  if(label!=NULL) free(label);
  if(vCatalog!=NULL) free(vCatalog);
  if(ux!=NULL) free(ux);
  if(uy!=NULL) free(uy);
  if(uxxx!=NULL) free(uxxx);
  if(uyyy!=NULL) free(uyyy);
  if(uxxy!=NULL) free(uxxy);
  if(uxyy!=NULL) free(uxyy);
  if(parVortex!=0) free(parVortex);
  if(Xload!=NULL) free(Xload);
  if(Yload!=NULL) free(Yload);
  if(Zload!=NULL) free(Zload);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

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