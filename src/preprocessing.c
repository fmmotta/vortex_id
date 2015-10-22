#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "preprocessing.h"

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int loadAxis(FILE *nFile,int Nx,int Ny,int Nz,
             double *X,double *Y,double *Z)
{
  const int Npre=20;
  int i,j,k,Nn;
  double a,b,c;
  char buffer[1024];

  if(nFile==NULL)
    return -3;
  
  for(i=0;i<Npre;i++)
    fgets(buffer,1024,nFile);
  fscanf(nFile,"%d",&Nn);
  
  fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
  X[0]=a;
  Y[0]=b;
  Z[0]=c;
  for(i=1;i<Nx+1;i+=1){
    fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    X[i]=a;
  }
  for(j=1;j<Ny+1;j+=1){
    fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    Y[j]=b;
    for(i=1;i<Nx+1;i+=1)
      fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
  }
  for(k=1;k<Nz+1;k+=1){
    fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    Z[k]=c;
    for(i=1;i<Nx+1;i+=1)
      fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    
    for(j=1;j<Ny+1;j+=1){
      fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
      for(i=1;i<Nx+1;i+=1)
        fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    }
  }

  return 0;
}

int printAxis(FILE *nFile,int Nx,int Ny,int Nz,
              double  *X,double  *Y,double *Z,
              double *X2,double *Y2,double *Z2,
              const char *folder)
{
  const int Npre=20;
  int i,j,k,err;
  char filename[100+1];
  FILE *ouFile;

  err=loadAxis(nFile,Nx,Ny,Nz,X,Y,Z);
  if(err!=0)
    return err;

  sprintf(filename,"%s/xAxis.dat",folder);
  ouFile = fopen(filename,"w");
  for(i=0;i<Nx;i+=1)
    fprintf(ouFile,"%.8g\n",(X[i]+X[i+1])/2.);
  fclose(ouFile);

  sprintf(filename,"%s/yAxis.dat",folder);
  ouFile = fopen(filename,"w");
  for(j=0;j<Ny;j+=1)
    fprintf(ouFile,"%.8g\n",(Y[j]+Y[j+1])/2.);
  fclose(ouFile);

  sprintf(filename,"%s/zAxis.dat",folder);
  ouFile = fopen(filename,"w");
  for(k=0;k<Nz;k+=1)
    fprintf(ouFile,"%.8g\n",(Z[k]+Z[k+1])/2.);
  fclose(ouFile);

  err=loadAxis(nFile,Nx,Ny-1,Nz,X2,Y2,Z2);
  if(err!=0)
    return err;

  sprintf(filename,"%s/y2Axis.dat",folder);
  ouFile = fopen(filename,"w");
  fprintf(ouFile,"%.8g\n",(Y2[0]+Y[Ny])/2.);
  for(j=1;j<Ny-1;j+=1)
    fprintf(ouFile,"%.8g\n",(Y2[j]+Y2[j+1])/2.);
  fclose(ouFile);

  return 0;
}

int loadFields(int Nx,int Ny,int Nz,FILE *uFile,FILE *pFile,
               openFoamIcoData *node)
{
  const int Npre=20;
  char buffer[1024];
  int i,j,k,Nu,Np;

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,uFile);
  fscanf(uFile,"%d",&Nu);

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,pFile);
  fscanf(pFile,"%d",&Np);
  
  if(Nu != Np || Nu != 2*(Nx*Ny*Nz)){
    fclose(uFile); fclose(pFile);
    return -10;
  }

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  // Nx*Ny*k+Nx*j+i
  /*
  for(i=0;i<Nx;i+=1){
    for(k=0;k<Nz;k+=1){
      for(j=0;j<Ny;j+=1){
        fscanf(pFile,"%lf",&(node[id(i,j,k)].p));
        fscanf(uFile," (%lf%lf%lf)",&(node[id(i,j,k)].u),
                                    &(node[id(i,j,k)].v),
                                    &(node[id(i,j,k)].w));
      }
    }
  }*/
  
  for(k=0;k<Nz;k+=1){
    for(j=0;j<Ny;j+=1){
      for(i=0;i<Nx;i+=1){
        fscanf(pFile,"%lf",&(node[id(i,j,k)].p));
        fscanf(uFile," (%lf%lf%lf)",&(node[id(i,j,k)].u),
                                    &(node[id(i,j,k)].v),
                                    &(node[id(i,j,k)].w));
      }
    }
  }

  fclose(uFile); fclose(pFile);

  return 0;
}

int printYZsplitPlanes(int Nx,int Ny, int Nz,openFoamIcoData *node,
                       double *X,double *Y,double *Z,const char *folder)
{
  int i,j,k;
  char filename[100+1];
  FILE *zFile,*yFile,*vFile;

  for(i=0;i<Nx;i+=1){
    sprintf(filename,"%s/z-%d.dat",folder,i);
    zFile=fopen(filename,"w");
    sprintf(filename,"%s/y-%d.dat",folder,i);
    yFile=fopen(filename,"w");
    sprintf(filename,"%s/plane-%d.dat",folder,i);
    vFile=fopen(filename,"w");
    
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1){
        fprintf(zFile,"%.8g\n",(Z[k]+Z[k+1])/2.0);
        fprintf(yFile,"%.8g\n",(Y[j]+Y[j+1])/2.0);
        fprintf(vFile,"%.8g %.8g\n",node[id(i,j,k)].w,
                                    node[id(i,j,k)].v);
      }

    fclose(zFile);
    fclose(yFile);
    fclose(vFile);
  }

  return 0;
}

int printXYsplitPlanes(int Nx,int Ny, int Nz,openFoamIcoData *node,
                       double *X,double *Y,double *Z,const char *folder)
{
  int i,j,k;
  char filename[100+1];
  FILE *zFile,*yFile,*vFile;

  printf("%.8g %.8g %.8g\n",node[id(0,0,0)].u,
                            node[id(0,0,0)].v,
                            node[id(0,0,0)].w);

  printf("%.8g %.8g %.8g\n",node[id(0,1,0)].u,
                            node[id(0,1,0)].v,
                            node[id(0,1,0)].w);

  printf("%.8g %.8g %.8g\n",node[id(1,0,0)].u,
                            node[id(1,0,0)].v,
                            node[id(1,0,0)].w);

  printf("%.8g %.8g %.8g\n",node[id(0,0,1)].u,
                            node[id(0,0,1)].v,
                            node[id(0,0,1)].w);
  for(k=0;k<Nz;k+=1){
    sprintf(filename,"%s/x-%d.dat",folder,k);
    zFile=fopen(filename,"w");
    sprintf(filename,"%s/y-%d.dat",folder,k);
    yFile=fopen(filename,"w");
    sprintf(filename,"%s/plane-%d.dat",folder,k);
    vFile=fopen(filename,"w");
    
    // #define id(i,j,k) (Nx*Ny*(k)+Nx*(j)+(i))
    for(j=0;j<Ny;j+=1)
      for(i=0;i<Nz;i+=1){
        fprintf(zFile,"%.8g\n",(X[i]+X[i+1])/2.0);
        fprintf(yFile,"%.8g\n",(Y[j]+Y[j+1])/2.0);
        fprintf(vFile,"%.8g %.8g\n",node[id(i,j,k)].u,
                                    node[id(i,j,k)].v);
      }

    fclose(zFile);
    fclose(yFile);
    fclose(vFile);
  }

  return 0;
}

int printYZcoordinates(int Nx,int Ny,int Nz,double *X,
                       double *Y,double *Z,const char *folder)
{
  int i,j,k;
  char filename[100+1];
  FILE *ouFile;

  sprintf(filename,"%s/coordinates.dat",folder);
  ouFile = fopen(filename,"w");
  i=0;
  for(k=0;k<Nz;k+=1)
    for(j=0;j<Ny;j+=1)
      fprintf(ouFile,"%.12f %.12f\n",(Z[k]+Z[k+1])/2.0,(Y[j]+Y[j+1])/2.0);
  fclose(ouFile);

  return 0;
}