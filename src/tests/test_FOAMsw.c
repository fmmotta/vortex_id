#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "preprocessing.h"
#include "floodFill.h"
#include "lambdaInit.h"
#include "stencilExtended.h"
#include "vortexExtraction.h"


#define fieldAlloc(ptr,size,type) ptr=(type*)malloc(size*sizeof(type)); \
                                  if(ptr==NULL){                        \
                                    printf("memory not allocked\n");    \
                                    return 1;                           \
                                  }                                     \
                                  else{                                 \
                                    for(i=0;i<size;i+=1)                \
                                      ptr[i]=(type) 0;                  \
                                  }                                     \

#define dbgPrint(i) printf("Hello - debug - %lf\n",(float)i)

int isThere(int i){
  if(i>0)
  	return 1;
  else
  	return 0;
}

int main(int argc,char** argv){
  const int padWidth=2, Pop=10;
  long int N=6;
  int Height, Width,Depth;
  int i,j,k,l,Npre,Nu,Np,Nx,Ny,Nz,Nn,err,auHeight,auWidth;
  int nbList[8],eqList[Pop],**eqClass,*label;//,*label;
  char buffer[1024],filename[100];
  double *sField=NULL,x,y,x0[2],dx[2];
  double Z[1000],X2[1000],Y2[1000],Z2[1000];
  double *gField=NULL,*g2Field=NULL,*uField=NULL,X[1000],Y[1000];
  double *uBuff=NULL,Xbuff[1000+4],Ybuff[1000+4];
  double *ux,*uy,*uxxy,*uxyy,*uxxx,*uyyy,*w,ua,ub;
  FILE *uFile,*pFile,*nFile,*ouFile;
  FILE *zFile,*yFile,*wFile,*vFile,*xFile,*uRFile;
  openFoamIcoData v[N],*node=NULL;

  if(argc!=4 && argc!=1){
    printf("wrong number of arguments, need exactly 2 input files, velocity and pressure files\n");
    return 1;
  }

  /*
  Nx=Depth;
  Nz=Width;
  Ny=Height;*/

  Ny = 96;
  Nx = 256;
  Nz = 192;

  Height = Ny;
  Width  = Nx;//Nz;
  Depth  = Nz;//Nx;
  
  dbgPrint(0);

  fieldAlloc(label ,Height*Width,int);
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
  
  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  dbgPrint(1);

  N = Nx*Ny*Nz;
  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  Npre=20; // preamble size

  if(argc==4){
    uFile = fopen(argv[1],"r"); // velocity file
    pFile = fopen(argv[2],"r"); // pressure file
    nFile = fopen(argv[3],"r"); // nodes positions file
  }
  else if(argc==1){
    uFile = fopen("data/DNS_OPEN_FOAM/22.1395/U","r");
    pFile = fopen("data/DNS_OPEN_FOAM/22.1395/p","r");
    nFile = fopen("data/DNS_OPEN_FOAM/constant/polyMesh/points","r");
  }

  if(uFile==NULL || pFile==NULL || nFile==NULL){
    printf("problems opening the files\n"); 
    return 1;
  }
  
  dbgPrint(2);

  err=loadAxis(nFile,Nx,Ny,Nz,X,Y,Z);
  if(err!=0)
    return err;
  fclose(nFile);

  /*
  xFile=fopen("data/x.txt");
  yFile=fopen("data/y.txt");
  uRFile=fopen("data/t22.1395_z064.dat","r");

  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1){
      fscanf(uRFile,"%lf%lf",&ua,&ub);
      uField[2*(i*Width+j)+0] = ua;
      uField[2*(i*Width+j)+1] = ub;
    }
  }*/

  dbgPrint(2.3);

  
  err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
  if(err!=0)
    printf("Problems with loadFields\n");

  dbgPrint(2.6);

  ouFile = fopen("data/mathematicaRefU.dat","w");
  k=64;
  for(j=0;j<Height;j+=1)
    for(i=0;i<Width;i+=1){
      uField[2*(j*Width+i)+0] = node[id(i,j,k)].v;
      uField[2*(j*Width+i)+1] = node[id(i,j,k)].u;
      fprintf(ouFile,"%lf %lf\n",node[id(i,j,k)].u,node[id(i,j,k)].v);
    }
  fclose(ouFile);ouFile=NULL;

  dbgPrint(3);

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");
  
  dbgPrint(4);

  err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
    printf("Problems in uFieldTouBuff\n");

  dbgPrint(5);

  // \partial_x \vec u
  err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
  if(err!=0)
    printf("Problems in UtoUx5point\n");

  dbgPrint(6);

  // \partial_y \vec u
  err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
  if(err!=0)
    printf("Problems in UtoUy5point\n");

  dbgPrint(7);
  

  dbgPrint(11);

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gField[4*(i*Width+j)+0] = ux[2*(i*Width+j)+0];
      gField[4*(i*Width+j)+1] = uy[2*(i*Width+j)+0];
      gField[4*(i*Width+j)+2] = ux[2*(i*Width+j)+1];
      gField[4*(i*Width+j)+3] = uy[2*(i*Width+j)+1];
    }
  
  dbgPrint(12);
  
  err = gradUtoLamb(Height,Width,gField,&sField);
  if(err!=0)
    printf("Problems in gradU2UtoLambda\n");
  
  dbgPrint(13);

  err = floodFill(sField,Width,Height,eqClass,label);
  if(err!=0)
    printf("Problems in floodFill\n");

  dbgPrint(14);

  err = renameLabels(Height,Width,label);
  if(err>0)
    printf("%d connected component(s)\n",err);
  else
    printf("problems with renameLabels - %d\n",err);

  dbgPrint(15);
  
  {
    FILE *dadosout;
    ouFile  =fopen("data/dens_z064.dat","w");
    dadosout=fopen("data/initFOAMsw.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        y = Y[i];
        x = X[j];
        
        fprintf(dadosout,"%f %f %.12f \n",x,y,sField[i*Width+j]);
        fprintf(ouFile,"%.12f\n",sField[i*Width+j]);
      }

    fclose(dadosout);dadosout=NULL;
    fclose(ouFile);ouFile=NULL;

    dadosout=fopen("data/labelFOAMsw.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = Y[i];
        x = X[j];
        
        fprintf(dadosout,"%f %f %2d \n",x,y,label[i*Width+j]+1);
      }
      fprintf(dadosout,"\n");
    }

    fclose(dadosout);

    dadosout=fopen("data/presentFOAMsw.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = Y[i];
        x = X[j];
        
        fprintf(dadosout,"%f %f %2d \n",x,y,isThere(label[i*Width+j]+1));
      }
      fprintf(dadosout,"\n");
    }

    fclose(dadosout);
  }
  
  dbgPrint(16);

  if(sField!=NULL)
    free(sField);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  return 0;
}
