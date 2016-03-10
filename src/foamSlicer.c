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
#include "mt64.h"
#include "ini.h"
#include "stencilExtended.h"
#include "lambdaInit.h"
#include "vortexGen.h"
#include "preprocessing.h"
#include "inputManager.h"
#include "vortexExtraction.h"		
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
                                  }           

int main(int argc,char **argv){
  int Height=100,Width=100,Depth;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum,Nsnapshots;
  int i,j,k,err,pln[8128],pn,n;
  double *uField,*Xload,*Yload,*Zload,*X,*Y,*Z,t0,dt,t;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1];
  FILE *dadosout,*uFile,*pFile,*nFile;
  openFoamIcoData *node=NULL;
  configVar cfg;
  
  planeNum=3;
  pln[0]=0; pln[1]=64; pln[2]=128;

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
    planeNum = 1;
    pln[0] = planeIndex;
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
  
  if(planeIndex<0)
    printf("Switching to plane number list\n");

  if(planeIndex>=Depth)
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
  
  dbgPrint(2,0);

  /**********************************/  

  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */
  
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(X,Nx,double);
  fieldAlloc(Y,Ny,double);
  fieldAlloc(Z,Nz,double);  
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
  if(err!=0){
  	printf("Coulnd't execute loadAxis\n");
    return err;
  }
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
    fprintf(dadosout,"%.8g\n",X[i]);
  fclose(dadosout);

  sprintf(filename,"%s/Yaxis.dat",folder);
  dadosout = fopen(filename,"w");
  for(i=0;i<Ny;i+=1)
    fprintf(dadosout,"%.8g\n",Y[i]);
  fclose(dadosout);

  sprintf(filename,"%s/Zaxis.dat",folder);
  dadosout = fopen(filename,"w");
  for(i=0;i<Nz;i+=1)
    fprintf(dadosout,"%.8g\n",Z[i]);
  fclose(dadosout);

  for(n=0;n<Nsnapshots;n+=1){
    
    printf("%d timesteps processed\n",n);
    
    t=t0+((double)n)*dt;
        
    for(i=0;i<2*Height*Width;i+=1)
      uField[i]=0.;
    
    dbgPrint(5,0);
    
    sprintf(filename,"%s/%g/U",foamFolder,t);
    uFile = fopen(filename,"r");
    if(uFile==NULL)
      printf("problems opening uFile - %s\n",filename);
    
    sprintf(filename,"%s/%g/p",foamFolder,t);
    pFile = fopen(filename,"r");
    if(pFile==NULL)
      printf("problems opening pFile - %s\n",filename);

    if(uFile==NULL || pFile==NULL)
      continue;
    
    dbgPrint(5,1);
    
    err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
    if(err!=0)
      printf("Problems with loadFields\n");
    
    fclose(pFile); fclose(uFile);
    
    dbgPrint(5,2);

    if(DEBUG_PRINT)
      printf("folder = %s\n",folder);
    
    dadosout=NULL;
    
    if(planeType==0){
      if(DEBUG_PRINT)
        printf("XY plane\n");
      
      if(planeIndex>=0){
        k=planeIndex;

      	sprintf(filename,"%s/plane-z%3d-%.4f.dat",folder,k,t);
      	dadosout=fopen(filename,"w");
        if(dadosout==NULL){
          printf("Could not open file\n");
          exit(EXIT_FAILURE);
        }

        for(j=0;j<Height;j+=1)
          for(i=0;i<Width;i+=1){
            uField[2*(j*Width+i)+0] = node[id(i,j,k)].u;
            uField[2*(j*Width+i)+1] = node[id(i,j,k)].v;
            fprintf(dadosout,"%.8g %.8g\n",uField[2*(j*Width+i)+0]
                                        ,uField[2*(j*Width+i)+1]);
          }
        
        if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
      }
      else{
      	for(pn=0;pn<planeNum;pn+=1){
      	  k=pln[pn];
      	  
      	  sprintf(filename,"%s/plane-z%3d-%.4f.dat",folder,k,t);
      	  dadosout=fopen(filename,"w");
          if(dadosout==NULL){
            printf("Could not open file\n");
            exit(EXIT_FAILURE);
          }

          for(j=0;j<Height;j+=1)
            for(i=0;i<Width;i+=1){
              uField[2*(j*Width+i)+0] = node[id(i,j,k)].u;
              uField[2*(j*Width+i)+1] = node[id(i,j,k)].v;
              fprintf(dadosout,"%.8g %.8g\n",uField[2*(j*Width+i)+0]
                                            ,uField[2*(j*Width+i)+1]);
            }
          
          if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
        }
      }
    }
    else if(planeType==1){
      if(DEBUG_PRINT)
        printf("YZ plane\n");
      
      if(planeIndex>=0){
        i=planeIndex;
        
        sprintf(filename,"%s/plane-x%3d-%.4f.dat",folder,planeIndex,t);
      	dadosout=fopen(filename,"w");
        if(dadosout==NULL){
          printf("Could not open file\n");
          exit(EXIT_FAILURE);
        }

        for(j=0;j<Height;j+=1)
          for(k=0;k<Width;k+=1){
            uField[2*(j*Width+k)+0] = node[id(i,j,k)].w;
            uField[2*(j*Width+k)+1] = node[id(i,j,k)].v;
            fprintf(dadosout,"%.8g %.8g\n",uField[2*(j*Width+k)+0]
                                          ,uField[2*(j*Width+k)+1]);
          }
        
        if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
      }
      else{
      	for(pn=0;pn<planeNum;pn+=1){
          i=pln[pn];
      	  
      	  sprintf(filename,"%s/plane-x%3d-%.4f.dat",folder,i,t);
      	  dadosout=fopen(filename,"w");
          if(dadosout==NULL){
            printf("Could not open file\n");
            exit(EXIT_FAILURE);
          }

          for(j=0;j<Height;j+=1)
            for(k=0;k<Width;k+=1){
              uField[2*(j*Width+k)+0] = node[id(i,j,k)].w;
              uField[2*(j*Width+k)+1] = node[id(i,j,k)].v;
              fprintf(dadosout,"%.8g %.8g\n",uField[2*(j*Width+k)+0]
                                            ,uField[2*(j*Width+k)+1]);
            }
          
          if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
        }
      }
    }
    else if(planeType==2){
      if(DEBUG_PRINT)
        printf("XZ plane\n");
      
      
      if(planeIndex>=0){
        j=planeIndex;
     
        sprintf(filename,"%s/plane-y%3d-%.4f.dat",folder,planeIndex,t);
        dadosout=fopen(filename,"w");
        if(dadosout==NULL){
          printf("Could not open file\n");
          exit(EXIT_FAILURE);
        }

        for(k=0;k<Height;k+=1)
          for(i=0;i<Width;i+=1){
            uField[2*(k*Width+i)+0] = node[id(i,j,k)].w;
            uField[2*(k*Width+i)+1] = node[id(i,j,k)].u;
            fprintf(dadosout,"%.8g %.8g\n",uField[2*(k*Width+i)+0]
                                          ,uField[2*(k*Width+i)+1]);
          }
        
        if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
      }
      else{
      	for(pn=0;pn<planeNum;pn+=1){
          j=pln[pn];
      	  
      	  sprintf(filename,"%s/plane-y%3d-%.4f.dat",folder,j,t);
      	  dadosout=fopen(filename,"w");
          if(dadosout==NULL){
            printf("Could not open file\n");
            exit(EXIT_FAILURE);
          }

          for(k=0;k<Height;k+=1)
            for(i=0;i<Width;i+=1){
              uField[2*(k*Width+i)+0] = node[id(i,j,k)].w;
              uField[2*(k*Width+i)+1] = node[id(i,j,k)].u;
              fprintf(dadosout,"%.8g %.8g\n",uField[2*(k*Width+i)+0]
                                            ,uField[2*(k*Width+i)+1]);
            }
          
          if(dadosout!=NULL){fclose(dadosout);dadosout=NULL;}
        }
      }
    }
    
    dbgPrint(5,3);
  }

  dbgPrint(6,1);

  free(uField);
  free(X);
  free(Y);
  free(Z);
  free(Xload);
  free(Yload);
  free(Zload);
  free(node);

  return 0;
}