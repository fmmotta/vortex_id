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

int vortexExtraction(int Height,int Width, int nCnect,
	                 float *x0, float *dx,float *sField,
	                 float *gField,int *label,float **vCatalogOut)
{
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2]; // vorticity
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float *vCatalog=NULL;

  if((Height<=0)||(Width<=0))
  	return -1;
  if(vCatalogOut==NULL)
    return -2;
  
  vCatalog = *vCatalogOut;
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
   *
   * update: stack overflow suggestion does not seems to work on my computer
   */

  /*   
   * BAD IDEA -- trocar essa bosta assim que possível
   * está a caminho
   *
  if(*vCatalogOut==NULL){
    vCatalog=(float*)malloc(4*nCnect*sizeof(float));
    if(vCatalog==NULL)
      return -3;
  }
  else{
    vCatalog=(float*)realloc(vCatalog,4*nCnect*sizeof(float));
    if(vCatalog==NULL)
      return -2;
  }
  */

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

        // I think that here is the factor of 1/2 
        // that was missing from my calculation
        // correcting and anotating for next commit
        // this shouldn't affect anything in the calculation
        // changing accordingly on the pdf
        A[k] += dx[0]*dx[1];
        w[k] += ( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)/1.12091; // Constant comming from lamb-oseen vortex;
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    G = 1.397948086*w[k]; // Constant comming from lamb-oseen vortex;

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}

int vortexExtSimple(int Height,int Width,float *x0, float *dx,
                    int **eqClass,float *sField,float *gField,int *label,
                    float threshold,int *nCnectOut,float **vCatalogOut){
  
  int i=0,err=0,pass=0,nCnect=0,nCnect0=0,it=0;
  float *vCatalog=NULL,majorVortex[4];
  vCatalog = *vCatalogOut;
  
  for(i=0;i<Height*Width;i+=1)
    label[i]=-1;

  err = gradUtoLamb(Height,Width,gField,&sField);
  if(err!=0)
    return err;
 
  err = floodFill(sField,Width,Height,eqClass,label);
  if(err!=0)
    return err;

  err = renameLabels(Height,Width,label);
  if(err>0)
    nCnect=err;
  else
    return err;

  err=vortexExtraction(Height,Width,nCnect,x0,dx,sField,
                       gField,label,&vCatalog);
  if(err!=0)
    return err; 

  *nCnectOut = nCnect;
  *vCatalogOut = vCatalog;

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

int greaterAbsCirculation(const float *v,const float *p){
  if(fabs(v[0])>fabs(p[0]))
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

int greaterAbsVorticity(const float *v,const float *p){
  if( fabs(v[0]/(v[1]*v[1]))>fabs(p[0]/(p[1]*p[1])) )
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

int vortexExtRecursive(int Height,int Width,float *x0, float *dx,int **eqClass,
                       float *sField,float *gField,int *label, float threshold, 
                       float *vCatalog, int *rCnectOut,float **rCatalogOut){
  int maxIt;
  int i=0,err=0,pass=0,rCnect=0,nCnect=0,nCnect0=0,it=0;
  float *rCatalog=NULL,majorVortex[4];

  if(Height<=0 || Width<=0)
    return -1;

  rCatalog = *rCatalogOut;

  do{
    // if(it>=maxIt) break;
    printf("it=%d\n",it);
    for(i=0;i<Height*Width;i+=1)
      label[i]=-1;

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0)
      return err;
  
    err = floodFill(sField,Width,Height,eqClass,label);
    if(err!=0)
      return err;

    err = renameLabels(Height,Width,label);
    if(err>0){
      nCnect=err;
    }
    else
      return err;

    err=vortexExtraction(Height,Width,nCnect,x0,dx,sField,
                         gField,label,&vCatalog);
    if(err!=0)
      return err;

    vortexQuickSort(vCatalog,nCnect,&greaterAbsCirculation);

    if(abs(vCatalog[4*0+0])>threshold){
      pass=1; 
      rCnect+=1;
      majorVortex[0]=-vCatalog[0]; rCatalog[4*it+0] = vCatalog[0];
      majorVortex[1]= vCatalog[1]; rCatalog[4*it+1] = vCatalog[1];
      majorVortex[2]= vCatalog[2]; rCatalog[4*it+2] = vCatalog[2];
      majorVortex[3]= vCatalog[3]; rCatalog[4*it+3] = vCatalog[3];
    }
    else
      break;

    err = addSingleOseen(1,majorVortex,x0,dx,Height,Width,&gField);
    if(err!=0){
      printf("alguma merda séria tá acontecendo\n");
      return err;}

    it+=1;
  }while(pass!=0);

  *rCnectOut = rCnect;
  *rCatalogOut = rCatalog; 

  return 0;
}

int applySwirlingStrengthThreshold(int Height,int Width,float *sField,
                                   float theta)
{
  int i,j;

  if(Height<=0 || Width<=0)
    return -1;

  if(sField==NULL)
    return -2;

  if(theta<0.)
    return -3;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      if(sField[i*Width+j]<= theta)
        sField[i*Width+j]=0.;
    }
  
  return 0;
}

double radiusRoot(double x, void *par)
{
  double val,*cut;
  val = (2.*x*x)/(exp(x*x)-1.0)-1.0; 
  cut = (double*)par;
  
  return (val-*cut);
}

int vortexExtTreshold(int Height,int Width, int nCnect,float theta,
                      float *x0, float *dx,float *sField,
                      float *gField,int *label,float **vCatalogOut)
{
  int i,j,k;
  float G,a,b,x,y,rc,rbar,gradU[2][2],eta; // vorticity
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float *vCatalog=NULL;
  int iter = 0, max_iter = 100,status;
  double x_lo = 0.1, x_hi = 1.1,cut=0.0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  
  if((Height<=0)||(Width<=0))
    return -1;
  if(vCatalogOut==NULL)
    return -2;
  
  vCatalog = *vCatalogOut;
  
  F.function = &radiusRoot;
    
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  //gsl_root_fsolver_set (s, &F, x_lo, x_hi);

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
        w[k] += ( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rbar = sqrt(A[k]/M_PI);
    a = a0[k]/w[k];
    b = b0[k]/w[k]; 
    w[k] = w[k]/A[k]; // average vorticity

    rbar = sqrt(A[k]/M_PI); // now find rc

    /*
     * use GSL root finding to solve
     * (2.*eta^2)/(exp(eta^2)-1.) = 1+(theta/w[k])^2
     */
    
    cut = (double) ( (theta*theta)/(w[k]*w[k]) );
    
    x_lo = 0.005, x_hi = 2.1;
    F.params = (void*)&cut;
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    do{
      iter++;
      status = gsl_root_fsolver_iterate (s);
      eta = (float) gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS)
        break;
    }while (status == GSL_CONTINUE && iter < max_iter);

    rc = rbar/eta; 

    G = (w[k]*A[k])/(1.-exp(-eta*eta));
    
    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b; 
  }

  gsl_root_fsolver_free (s); // probably make a clause as to
                             // free this solver 
  *vCatalogOut=vCatalog;

  return 0;
}

int vortexExt2ndSwirl(int Height,int Width, int nCnect,
                      float *x0, float *dx,float *sField,
                      float *gField,int *label,float **vCatalogOut)
{
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2]; // vorticity
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
        
        x = x0[0] + i*dx[0];
        y = x0[1] + j*dx[1];

        A[k] += dx[0]*dx[1];
        w[k] += ( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrtf(A[k]/M_PI)*sqrtf(2.); // Constant comming from lamb-oseen vortex;
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    G = 0.8243606353500641*rc*rc*w[k]; // Constant comming from lamb-oseen 
                                       // vortex; Based on laplacian of omega
                                       // equals to sqrt(e)/2
    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}