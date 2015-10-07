#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

#define DEBUG_MODE false

#define id(i,j,k) (Nx*Ny*(k)+Nx*(j)+(i))

typedef struct openFoamIcoData {
    float u,v,w,p;
} openFoamIcoData;

int comp (const void * elem1, const void * elem2);

int main(int argc,char** argv){
  long int N=6;
  int i,j,k,l,Npre,Nu,Np,Nx,Ny,Nz,Nn;
  char buffer[1024];
  char filename[100];
  float dx,dz,Y[1000],dU,dV,dW;
  FILE *uFile,*pFile,*nFile,*ouFile;
  FILE *zFile,*yFile,*wFile,*vFile;
  openFoamIcoData v[N],*node=NULL;

  if(argc!=4){
    printf("wrong number of arguments, need exactly 2 input files, velocity and pressure files");
    return 1;
  }

  Nx=256;//Nx=256;
  Nz=192;//Nz=192;
  Ny=96;//Ny=96;

  dx = 2.0*M_PI/((float)Nx);
  dz = 1.0*M_PI/((float)Nz);

  Npre=20; // preamble size

  uFile = fopen(argv[1],"r"); // velocity file
  pFile = fopen(argv[2],"r"); // pressure file
  nFile = fopen(argv[3],"r"); // nodes positions file

  if(uFile==NULL || pFile==NULL || nFile==NULL){
    printf("problems opening the files\n"); 
    return 1;
  }

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,uFile);
  fscanf(uFile,"%d",&Nu);

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,pFile);
  fscanf(pFile,"%d",&Np);
  
  for(i=0;i<Npre;i++)
    fgets(buffer,1024,nFile);
  fscanf(nFile,"%d",&Nn);

  if(Nu != Np || Nu != 2*(Nx*Ny*Nz)){
  	printf("Non-Matching number of elements - Probably wrong pair of files\n"
           "Or wrong number of cell sizes\n"); 
  	fclose(uFile); fclose(pFile); fclose(nFile);
    return 2;
  }

  fgets(buffer,1024,nFile);
  fgets(buffer,1024,nFile);

  for(j=0;j<Ny+1;j+=1){
    float a,b,c;
    fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    Y[j]=b;
    for(i=0;i<Nx+1;i+=1)
      fscanf(nFile," (%lf%lf%lf)",&a,&b,&c);
    
    if(DEBUG_MODE)
      printf("Y[%d]=%lf\n",j,Y[j]);
  }
  fclose(nFile);
  //qsort(Y,Ny+1,sizeof(int),comp); // qsort is introducing some kind of bug

  ouFile = fopen("yAxis.dat","w");
  for(j=0;j<Ny;j+=1)
    fprintf(ouFile,"%lf\n",(Y[j]+Y[j+1])/2.);
  fclose(ouFile);

  ouFile = fopen("xAxis.dat","w");
  for(k=0;k<Nz;k+=1)
    fprintf(ouFile,"%lf\n",(k+0.5)*dz);
  fclose(ouFile);

  N = Nx*Ny*Nz;
  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  for(i=0;i<Nx;i+=1){
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1)
        fscanf(uFile," (%lf%lf%lf)",&(node[id(i,j,k)].u),
                                    &(node[id(i,j,k)].v),
                                    &(node[id(i,j,k)].w));
    
  }

  if(DEBUG_MODE){
    ouFile = fopen("check.dat","w");
    for(i=0;i<Nx;i+=1)
      for(k=0;k<Nz;k+=1)
        for(j=0;j<Ny;j+=1)
          fprintf(ouFile,"%f %f %f %f %f %f\n",(i+0.5)*dx,(k+0.5)*dz,
                                               (Y[j]+Y[j+1])/2.0,
                                               node[id(i,j,k)].u,
                                               node[id(i,j,k)].v,
                                               node[id(i,j,k)].w);
    fclose(ouFile);
  }
  
  for(i=0;i<Nx;i+=1){
    sprintf(filename,"planes/z-%d.dat",i);
    zFile=fopen(filename,"w");
    sprintf(filename,"planes/y-%d.dat",i);
    yFile=fopen(filename,"w");
    sprintf(filename,"planes/w-%d.dat",i);
    wFile=fopen(filename,"w");
    sprintf(filename,"planes/v-%d.dat",i);
    vFile=fopen(filename,"w");

    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1){
        fprintf(zFile,"%.8g\n",(k+0.5)*dz);
        fprintf(yFile,"%.8g\n",(Y[j]+Y[j+1])/2.0);
        fprintf(wFile,"%.8g\n",node[id(i,j,k)].w);
        fprintf(vFile,"%.8g\n",node[id(i,j,k)].v);
      }
      /*  fprintf(ouFile,"%.8g %.8g %.8g %.8g\n",(k+0.5)*dz,
                                               (Y[j]+Y[j+1])/2.0,
                                               node[id(i,j,k)].w,
                                               node[id(i,j,k)].v);*/
    fclose(ouFile);  
  }
  
  ouFile = fopen("coordinates.dat","w");
  i=0;
  for(k=0;k<Nz;k+=1)
    for(j=0;j<Ny;j+=1){
      if(DEBUG_MODE){
        printf("k=%d j=%d\n",k,j);
        printf("z=%f y=%f\n",(k+0.5)*dz,(Y[j]+Y[j+1])/2.0);
      }
      fprintf(ouFile,"%.12f %.12f\n",(k+0.5)*dz,
                                     (Y[j]+Y[j+1])/2.0);
  }
  fclose(ouFile);

  fclose(uFile); fclose(pFile);

  return 0;
}

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

/*

  sprintf(filename,"planes/zCheck.dat");
  ouFile=fopen(filename,"w");  
  for(i=0;i<Nx;i+=1){
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1)
        fprintf(ouFile,"%.8g %.8g %.8g\n",//(k+0.5)*dz,
                                          //(Y[j]+Y[j+1])/2.0,
                                      node[id(i,j,k)].u,
                                      node[id(i,j,k)].v,
                                      node[id(i,j,k)].w);
    
  }
  fclose(ouFile);
*/
 /*

  for(i=0;i<Nx;i+=1){
    sprintf(filename,"planes/xCut-%d.dat",i);
    ouFile=fopen(filename,"w");
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Ny;j+=1)
        fprintf(ouFile,"%.8g %.8g\n",//(k+0.5)*dz,
                                     //(Y[j]+Y[j+1])/2.0,
                                      node[id(i,j,k)].w,
                                      node[id(i,j,k)].v);
    fclose(ouFile);  
  }

 */