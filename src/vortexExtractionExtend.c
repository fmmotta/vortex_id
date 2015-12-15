#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

inline int add_dgradU(int Height,int Width, int i,int j,int ik,int jk,
                      double *gField,double dgradU[][2],double X[],double Y[])
{
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      dgradU[0][1] = 9.*gField[4*(     i*Width +  j  )+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +  j  )+2];
    }
    else{
      dgradU[0][1] = 9.*gField[4*(     i*Width   +   j  )+1]+
                     3.*gField[4*(     i*Width   +(j+jk))+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*(     i*Width +(j+jk))+2];
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      dgradU[0][1] = 9.*gField[4*(     i*Width     +  j )+1]+
                     3.*gField[4*((i+ik)*Width     +  j )+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*((i+ik)*Width +   j  )+2];
    }
    else{
      dgradU[0][1] = 9.*gField[4*(     i*Width     +  j )+1]+
                     3.*gField[4*((i+ik)*Width     +  j )+1]+
                     3.*gField[4*(     i*Width   +(j+jk))+1]+
                     1.*gField[4*((i+ik)*Width   +(j+jk))+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*((i+ik)*Width +   j  )+2]+
                     3.*gField[4*(     i*Width +(j+jk))+2]+
                     1.*gField[4*((i+ik)*Width +(j+jk))+2];  
    }
  }

  dgradU[0][1] *= fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0;
  dgradU[1][0] *= fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0;

  return 0;
}

inline int add_dxdy(int Height,int Width, int i,int j,int ik,int jk,
                    double *gField,double *dx,double *dy,double X[],double Y[])
{
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *dy = 9.*gField[4*(     i*Width  + j  )+2]*Y[i];

      *dy-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i];
        
      *dx = 9.*gField[4*(     i*Width  + j  )+2]*X[j];

      *dx-= 9.*gField[4*(     i*Width  + j  )+1]*X[j];
    }
    else{
      *dy = 9.*gField[4*(     i*Width  + j  )+2]*Y[i]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i];

      *dy-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i];
        
      *dx = 9.*gField[4*(     i*Width  + j  )+2]*X[j]+
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk];

      *dx-= 9.*gField[4*(     i*Width  + j  )+1]*X[j]+
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk];
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      *dy = 9.*gField[4*(     i*Width  + j  )+2]*Y[i]   +
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik];

      *dy-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i]   +
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik];
        
      *dx = 9.*gField[4*(     i*Width  + j  )+2]*X[j]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j];

      *dx-= 9.*gField[4*(     i*Width  + j  )+1]*X[j]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j];
    }
    else{
      *dy = 9.*gField[4*(     i*Width  + j  )+2]*Y[i]   +
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i]   +
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*Y[i+ik];

      *dy-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i]   +
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i]   +
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*Y[i+ik];
        
      *dx = 9.*gField[4*(     i*Width  + j  )+2]*X[j]   +
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j]   +
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk]+
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*X[j+jk];

      *dx-= 9.*gField[4*(     i*Width  + j  )+1]*X[j]   +
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j]   +
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk]+
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*X[j+jk];
    }
  }

  *dy *= fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0;
  *dx *= fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0;

  return 0;
}

inline int add_dA(int Height,int Width, int i,int j,int k,int ik,int jk,
                  double *gField,int *label,double *dA,double X[],double Y[])
{
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0);
    }
    else{
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[     i*Width+(j+jk)]==k)?1.0:0.0);
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[(i+ik)*Width  +  j ]==k)?1.0:0.0);
    }
    else{
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[(i+ik)*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[     i*Width+(j+jk)]==k)?1.0:0.0)+
            1.*((label[(i+ik)*Width+(j+jk)]==k)?1.0:0.0);
    }
  }
    
       
  *dA *= fabs((Y[i+1]-Y[i])*(X[j+1]-X[j]))/64.0;

  return 0;
}

int vortexExtFromVortCurv(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut)
{
  int i,j,k,err;
  double G,a,b,rc,dx,dy,dA,dgradU[2][2]; // vorticity
  double w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  double *vCatalog=NULL,cutoff=0.001;

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
        // ++ quadrant        
        err=add_dgradU(Height,Width,i,j,1,1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    // 0.977816 corrects for gridsize
    rc= sqrt(A[k]/M_PI)*(sqrtf(2.));//0.977816); // Constant comming from 
                                                 //  lamb-oseen vortex;
    if(fabs(w[k])>0.){
      a=a0[k]/w[k]; 
      b=b0[k]/w[k];
    }
    else{
      a=X[0];
      b=Y[0];
    }
    
    G = 2.541494083*w[k]; // 2.541494083 = 1/(1-1/sqrt(e)) 
                          // ... should correct for finite grid size?
    
    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}

int vortexExtFromSwirlStr(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut)
{
  int i,j,k,err;
  double G,a,b,rc,dx,dy,dA,dgradU[2][2]; // vorticity
  double w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  double *vCatalog=NULL;

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
        // ++ quadrant        
        err=add_dgradU(Height,Width,i,j,1,1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        A[k] += dA;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    // 0.977816 corrects for gridsize
    rc= sqrt(A[k]/M_PI)/1.12091;//0.977816); // Constant comming from 
                                                 //  lamb-oseen vortex;
    if(fabs(w[k])>0.){
      a=a0[k]/w[k]; 
      b=b0[k]/w[k];
    }
    else{
      a=X[0]-1.;
      b=Y[0]-1.;
    }
        
    G = 1.397948086*w[k];
    
    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}