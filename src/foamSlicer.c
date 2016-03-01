#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "ini.h"
#include "preprocessing.h"
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
                                  }           

int main(int argc,char **argv){
  int Height=100,Width=100,Depth;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum;
  int i,j,k,err,pln[256];
  double *uField,*Xload,*Yload,*Zload,*X,*Y,*Z;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1];
  FILE *dadosin,*dadosout,*uFile,*pFile,*nFile;
  openFoamIcoData *node=NULL;
  
  planeNum=3;
  pln[0]=0; pln[1]=64; pln[2]=128;

  if(argc!=3){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  err=initConfig(&cfg);
  
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
  
  dbgPrint(0,0);

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

  dbgPrint(1,0);

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }
  
  runType = cfg.runType;
  
  dbgPrint(2,0);

  /**********************************/  

  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */
  
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(Xload,Nx+1,double);
  fieldAlloc(Yload,Ny+1,double);
  fieldAlloc(Zload,Nz+1,double);
  
  dbgPrint(3,0);

  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  dbgPrint(4,0);

  sprintf(filename,"%s/constant/polyMesh/points",foamFolder);
  nFile = fopen(filename,"r");
  err=loadAxis(nFile,Nx,Ny,Nz,Xload,Yload,Zload);
  if(err!=0)
    return err;
  fclose(nFile);
  
  for(j=0;j<Nx;j+=1)
    X[j] = (Xload[j]+Xload[j+1])/2.;
  for(i=0;i<Ny;i+=1)
    Y[i] = (Yload[i]+Yload[i+1])/2.;
  for(k=0;k<Nz;k+=1)
    Z[k] = (Zload[k]+Zload[k+1])/2.;
  
  sprintf(filename,"%s/Xaxis.dat",folder);
  dadosout = fopen(filename,"w");
  for(i=0;i<Nx;i+=1)
    fprintf(dadosout,"%f\n",X[i]);
  fclose(dadosout);

  sprintf(filename,"%s/Yaxis.dat",folder);
  dadosout = fopen(filename,"w");
  for(i=0;i<Ny;i+=1)
    fprintf(dadosout,"%f\n",Y[i]);
  fclose(dadosout);

  sprintf(filename,"%s/Zaxis.dat",folder);
  dadosout = fopen(filename,"w");
  for(i=0;i<Nz;i+=1)
    fprintf(dadosout,"%f\n",Z[i]);
  fclose(dadosout);

  for(n=0;n<Nsnapshots;n+=1){
    
    printf("%d runs have passed\n",n);

    t=t0+((double)n)*dt;
        
    for(i=0;i<2*Height*Width;i+=1)
      uField[i]=0.;
    
    dbgPrint(5,0);

    sprintf(filename,"%s/%.4f/U",foamFolder,t);
    uFile = fopen(filename,"r");
    if(uFile==NULL)
      printf("problems opening uFile - %d\n",n);

    sprintf(filename,"%s/%.4f/p",foamFolder,t);
    pFile = fopen(filename,"r");
    if(pFile==NULL)
      printf("problems opening pFile - %d\n",n);
    
    dbgPrint(5,1);

    err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
    if(err!=0)
      printf("Problems with loadFields\n");
    
    fclose(pFile); fclose(uFile);
    
    dbgPrint(5,2);

    if(DEBUG_PRINT)
      printf("folder = %s\n",folder);
    
    if(planeType==0)
      sprintf(filename,"%s/plane-z%3d-%.4f.dat",folder,pln[0],t);
    dadosout=NULL;
    dadosout=fopen(filename,"w");
    if(dadosout==NULL){
      printf("Could not open file\n");
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
          fprintf(dadosout,"%lf %lf\n",uField[2*(j*Width+i)+0]
                                      ,uField[2*(j*Width+i)+1]);
        }
    }
    else if(planeType==2){
      if(DEBUG_PRINT)
        printf("XZ plane\n");
      j=planeIndex;
      for(k=0;k<Height;k+=1)
        for(i=0;i<Width;i+=1){
          uField[2*(k*Width+i)+0] = node[id(i,j,k)].w;
          uField[2*(k*Width+i)+1] = node[id(i,j,k)].v;
          fprintf(dadosout,"%lf %lf\n",uField[2*(j*Width+i)+0]
                                      ,uField[2*(j*Width+i)+1]);
        }
    }
    
    if(dadosout!=NULL)
     fclose(dadosout);

    dbgPrint(5,3);
  }

  return 0;
}