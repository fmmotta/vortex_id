#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

int vortexExtraction(int Height,int Width, int nCnect,
	                 float *x0, float *dx,float *sField,
	                 float *gField,int *label,float **vCatalogOut)
{
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2]; // vorticity
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float *vCatalog;

  if((Height<=0)||(Width<=0))
  	return -1;

  if(vCatalogOut==NULL){
    vCatalog=(float*)malloc(4*nCnect*sizeof(float));
    if(vCatalog==NULL)
      return -3;
  }
  else{
    /*
     * WARNING : 
     *
     * potentially problematic behaviour, but works with
     * glibc, but have to check if it works for other libs.
     * I'm exploring the fact that if the size of the allocated array
     * is the same as the one I'm resizing it to, 
     * realloc changes nothing. 
     * see: http://stackoverflow.com/questions/18617620/
     *      behavior-of-realloc-when-the-new-size-is-the-same-as-the-old-one
     */
    vCatalog=(float*)realloc(*vCatalogOut,4*nCnect*sizeof(float));
    if(vCatalog==NULL)
      return -2;
  }

  for(k=0;k<nCnect;k+=1){
    w[k]=0.;A[k]=0.;a0[k]=0.;b0[k]=0.;
  }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){

        gradU[0][0] = gField[4*(i*Width+j)+0];
        gradU[0][1] = gField[4*(i*Width+j)+1];
        gradU[1][0] = gField[4*(i*Width+j)+2];
        gradU[1][1] = gField[4*(i*Width+j)+3];
        
        x = x0[0] + i*dx[0];
        y = x0[1] + j*dx[1];

        A[k] += dx[0]*dx[1];
        w[k] += ((gradU[1][0]-gradU[0][1])/2.)*dx[0]*dx[1];
        a0[k] += x*((gradU[1][0]-gradU[0][1])/2.)*dx[0]*dx[1];
        b0[k] += y*((gradU[1][0]-gradU[0][1])/2.)*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)/1.12091; // Constant comming from lamb-oseen vortex;
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    G = 2.7958961719283355*w[k]; // Constant comming from lamb-oseen vortex;

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}

/*
 * Quicksort implementation inspired on Roseta Code:
 * http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
 */

void vortexQuickSort(float *v,int nCnect,
                     int (*cmp)(const float*,const float*)){
  int i,j;
  float p[4],t[4];
  if(nCnect<2)
    return;

  p[0] = v[4*(nCnect/2)+0];
  p[1] = v[4*(nCnect/2)+1];
  p[2] = v[4*(nCnect/2)+2];
  p[3] = v[4*(nCnect/2)+3];
  for(i=0,j=nCnect-1;;i++,j--){
    //while( v[4*i+0]/(v[4*i+1]*v[4*i+1]) > p[0]/(p[1]*p[1]) )
    while( cmp(v+4*i,p) )
      i++;
    //while( p[0]/(p[1]*p[1]) > v[4*j+0]/(v[4*j+1]*v[4*j+1]) )
    while( cmp(p,v+4*j) )
      j--;
    if(i >= j)
      break;

    t[0]=v[4*i+0]; t[1]=v[4*i+1];
    t[2]=v[4*i+2]; t[3]=v[4*i+3];

    v[4*i+0]=v[4*j+0]; v[4*i+1]=v[4*j+1]; 
    v[4*i+2]=v[4*j+2]; v[4*i+3]=v[4*j+3];

    v[4*j+0]=t[0]; v[4*j+1]=t[1]; 
    v[4*j+2]=t[2]; v[4*j+3]=t[3];
  }
  vortexQuickSort(v,i,cmp);
  vortexQuickSort(v+4*i,nCnect-i,cmp);
}

int lesserCirculation(const float *v,const float *p){
  if(v[0]<p[0])
    return 1;
  else
    return 0;
}

int greaterCirculation(const float *v,const float *p){
  if(v[0]>p[0])
    return 1;
  else
    return 0;
}

int lesserVorticity(const float *v,const float *p){
  if(v[0]/(v[1]*v[1])<p[0]/(p[1]*p[1]))
    return 1;
  else
    return 0;
}

int greaterVorticity(const float *v,const float *p){
  if(v[0]/(v[1]*v[1])>p[0]/(p[1]*p[1]))
    return 1;
  else
    return 0;
}

int lesserRadius(const float *v,const float *p){
  if(v[1]<p[1])
    return 1;
  else
    return 0;
}

int greaterRadius(const float *v,const float *p){
  if(v[1]>p[1])
    return 1;
  else
    return 0;
}

int vortexExtRecursive(int Height,int Width, int nCnect,
                       float *x0, float *dx,int **eqClass,
                       float *sField,float *gField,int *label,
                       float *rCatalogOut){
  int i,err=0,pass=0;
  float *vCatalog;
  do{
    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0)
      printf("Problems in gradUtoLamb\n");
  
    err = floodFill(sField,Width,Height,eqClass,label);
    if(err!=0)
      printf("Problems in floodFill\n");

    err = renameLabels(Height,Width,label);
    if(err>0)
      nCnect=err;
    else
      printf("problems with renameLabels - %d\n",err);

    err=vortexExtraction(Height,Width,nCnect,x0,dx,sField,
                         gField,label,&vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }
    
    
  }while(pass!=0);

  return 0;
}