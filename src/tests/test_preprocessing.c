#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "preprocessing.h"

int main(int argc,char** argv){
  long int N=6;
  int i,j,k,l,Npre,Nu,Np,Nx,Ny,Nz,Nn,err;
  char buffer[1024],filename[100];
  double dx,dz,Y[1000],X[1000],Z[1000],dU,dV,dW;
  double X2[1000],Y2[1000],Z2[1000];
  FILE *uFile,*pFile,*nFile,*ouFile;
  FILE *zFile,*yFile,*wFile,*vFile;
  openFoamIcoData v[N],*node=NULL;

  if(argc!=4 && argc!=1){
    printf("wrong number of arguments, need exactly 2 input files, velocity and pressure files\n");
    return 1;
  }

  Nx=256;
  Nz=192;
  Ny=96; 

  dx = 2.0*M_PI/((double)Nx);
  dz = 1.0*M_PI/((double)Nz);

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

  err=printAxis(nFile,Nx,Ny,Nz,X,Y,Z,X2,Y2,Z2,"data");
  if(err!=0)
    printf("problems with printAxis\n");
  fclose(nFile);
  
  N = Nx*Ny*Nz;
  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
  if(err!=0)
    printf("Problems with loadFields\n");

  err=printXYsplitPlanes(Nx,Ny,Nz,node,X,Y,Z,"data/planes");
  if(err!=0)
    printf("Problems on printYZsplitPlanes\n");
  /*
  err=printYZcoordinates(Nx,Ny,Nz,X,Y,Z,"data/planes");
  if(err!=0)
    printf("printYZcoordinates\n");
  /*
  err=printUField(Nx,Ny,Nz,node,X,Y,Z,"data/planes");
  if(err!=0)
    printf("printYZcoordinates\n");*/

  return 0;
}