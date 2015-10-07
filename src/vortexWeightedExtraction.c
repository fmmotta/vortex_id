#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

/*
    Weight calculation for the uniform 2D integration stencil, which have 
    the form :

       Dx Dy  [1  6 1]
      ------- [6 36 6]
         64   [1  6 1]
*/

float weight_neighbours(int i,int j,int *label,int Width,int Height){
  float weight=36.0; // counting the weight of the point itseft

  // 4-way conectivity
  if(bound_check(i+1,j,Width,Height) && (label[(i+1)*Width+j]>=0))
    weight+=6.0;  

  if(bound_check(i,j+1,Width,Height) && (label[(i)*Width+(j+1)]>=0))
    weight+=6.0;
  
  if(bound_check(i-1,j,Width,Height) && (label[(i-1)*Width+(j)]>=0))
    weight+=6.0;
  
  if(bound_check(i,j-1,Width,Height) && (label[(i)*Width+(j-1)]>=0))
    weight+=6.0;
  
  // 8-way conectivity
  if(bound_check(i+1,j+1,Width,Height) && (label[(i+1)*Width+(j+1)]>=0))
    weight+=1.0;
  
  if(bound_check(i+1,j-1,Width,Height) && (label[(i+1)*Width+(j-1)]>=0))
    weight+=1.0;
  
  if(bound_check(i-1,j+1,Width,Height) && (label[(i-1)*Width+(j+1)]>=0))
    weight+=1.0;
  
  if(bound_check(i-1,j-1,Width,Height) && (label[(i-1)*Width+(j-1)]>=0))
    weight+=1.0;
  
  return weight/64.0;
}

/*
  This weighted extraction is certanly giving sub-optimal results
  comparing with the default 'dumb' calculation.

  Warning : Be cautious, maybe this routine is systematically 
  giving bad results
 */

int vortexUnifWeightedExtraction(int Height,int Width, int nCnect,
                                 float *x0, float *dx,float *sField,
                                 float *gField,int *label,float **vCatalogOut)
{ 
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2],weight; 
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float *vCatalog=NULL;

  if((Height<=0)||(Width<=0))
    return -1;
  if(vCatalogOut==NULL)
    return -2;

  vCatalog = *vCatalogOut;
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
        
        y = x0[0] + i*dx[0]; 
        x = x0[1] + j*dx[1]; 

        weight = weight_neighbours(i,j,label,Width,Height);

        A[k] += weight*dx[0]*dx[1];
        w[k] += weight*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += weight*x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += weight*y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)/1.12091; 
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    G = 1.397948086*w[k]; 

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}

/*

int vortexUnifWeightedExtraction(int Height,int Width, int nCnect,
                                 float *x0, float *dx,float *sField,
                                 float *gField,int *label,float **vCatalogOut)
{ 
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2],weight; 
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float *vCatalog=NULL;

  if((Height<=0)||(Width<=0))
    return -1;
  if(vCatalogOut==NULL)
    return -2;

  vCatalog = *vCatalogOut;
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
        
        y = x0[0] + i*dx[0]; 
        x = x0[1] + j*dx[1]; 

        weight = weight_neighbours(i,j,label,Width,Height);

        A[k] += weight*dx[0]*dx[1];
        w[k] += weight*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += weight*x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += weight*y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)/1.12091; 
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    G = 1.397948086*w[k]; 

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
} */