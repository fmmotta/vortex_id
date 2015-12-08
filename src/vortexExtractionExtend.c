#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

int vortexExtractionExtend(int Height,int Width, int nCnect,double *X,double *Y,
                           double *sField,double *gField,int *label,
                           double **vCatalogOut)
{
  int i,j,k;
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

  for(i=1;i<Height-1;i+=1)
    for(j=1;j<Width-1;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant

        dgradU[0][1] = 9.*gField[4*(    i*Width     +j  )+1]+
                       3.*gField[4*((i+1)*Width     +j  )+1]+
                       3.*gField[4*(    i*Width   +(j+1))+1]+
                       1.*gField[4*((i+1)*Width   +(j+1))+1];

        dgradU[1][0] = 9.*gField[4*(    i*Width     +j  )+2]+
                       3.*gField[4*((i+1)*Width     +j  )+2]+
                       3.*gField[4*(    i*Width   +(j+1))+2]+
                       1.*gField[4*((i+1)*Width   +(j+1))+2];

        dgradU[0][1] *= (Y[i+1]-Y[i])*(X[j+1]-X[j])/64.0;
        dgradU[1][0] *= (Y[i+1]-Y[i])*(X[j+1]-X[j])/64.0;

        w[k] += dgradU[1][0]-dgradU[0][1];
        
        dy = 9.*gField[4*(    i*Width  +j  )+2]*Y[i]  +
             3.*gField[4*((i+1)*Width  +j  )+2]*Y[i+1]+
             3.*gField[4*(    i*Width+(j+1))+2]*Y[i]  +
             1.*gField[4*((i+1)*Width+(j+1))+2]*Y[i+1];

        dy-= 9.*gField[4*(    i*Width  +j  )+1]*Y[i]  +
             3.*gField[4*((i+1)*Width  +j  )+1]*Y[i+1]+
             3.*gField[4*(    i*Width+(j+1))+1]*Y[i]  +
             1.*gField[4*((i+1)*Width+(j+1))+1]*Y[i+1];
        
        dx = 9.*gField[4*(    i*Width  +j  )+2]*X[j]  +
             3.*gField[4*((i+1)*Width  +j  )+2]*X[j]  +
             3.*gField[4*(    i*Width+(j+1))+2]*X[j+1]+
             1.*gField[4*((i+1)*Width+(j+1))+2]*X[j+1];

        dx-= 9.*gField[4*(    i*Width  +j  )+1]*X[j]  +
             3.*gField[4*((i+1)*Width  +j  )+1]*X[j]  +
             3.*gField[4*(    i*Width+(j+1))+1]*X[j+1]+
             1.*gField[4*((i+1)*Width+(j+1))+1]*X[j+1];
        
        dy *= (Y[i+1]-Y[i])*(X[j+1]-X[j])/64.0;
        dx *= (Y[i+1]-Y[i])*(X[j+1]-X[j])/64.0;
        
        a0[k] += dx;
        b0[k] += dy;

        dA = 9.*((label[    i*Width  +  j]==k)?1.0:0.0)+
             3.*((label[(i+1)*Width  +  j]==k)?1.0:0.0)+
             3.*((label[    i*Width+(j+1)]==k)?1.0:0.0)+
             1.*((label[(i+1)*Width+(j+1)]==k)?1.0:0.0);

        A[k] += dA*(Y[i+1]-Y[i])*(X[j+1]-X[j])/64.0;
        
        /*************************************************/
        
        // -+ quadrant

        dgradU[0][1] = 9.*gField[4*(    i*Width     +j  )+1]+
                       3.*gField[4*((i-1)*Width     +j  )+1]+
                       3.*gField[4*(    i*Width   +(j+1))+1]+
                       1.*gField[4*((i-1)*Width   +(j+1))+1];

        dgradU[1][0] = 9.*gField[4*(    i*Width     +j  )+2]+
                       3.*gField[4*((i-1)*Width     +j  )+2]+
                       3.*gField[4*(    i*Width   +(j+1))+2]+
                       1.*gField[4*((i-1)*Width   +(j+1))+2];

        dgradU[0][1] *= (Y[i]-Y[i-1])*(X[j+1]-X[j])/64.0;
        dgradU[1][0] *= (Y[i]-Y[i-1])*(X[j+1]-X[j])/64.0;

        w[k] += dgradU[1][0]-dgradU[0][1];
        
        dy = 9.*gField[4*(    i*Width  +j  )+2]*Y[i]  +
             3.*gField[4*((i-1)*Width  +j  )+2]*Y[i-1]+
             3.*gField[4*(    i*Width+(j+1))+2]*Y[i]  +
             1.*gField[4*((i-1)*Width+(j+1))+2]*Y[i-1];

        dy-= 9.*gField[4*(    i*Width  +j  )+1]*Y[i]  +
             3.*gField[4*((i-1)*Width  +j  )+1]*Y[i-1]+
             3.*gField[4*(    i*Width+(j+1))+1]*Y[i]  +
             1.*gField[4*((i-1)*Width+(j+1))+1]*Y[i-1];
        
        dx = 9.*gField[4*(    i*Width  +j  )+2]*X[j]  +
             3.*gField[4*((i-1)*Width  +j  )+2]*X[j]  +
             3.*gField[4*(    i*Width+(j+1))+2]*X[j+1]+
             1.*gField[4*((i-1)*Width+(j+1))+2]*X[j+1];

        dx-= 9.*gField[4*(    i*Width  +j  )+1]*X[j]  +
             3.*gField[4*((i-1)*Width  +j  )+1]*X[j]  +
             3.*gField[4*(    i*Width+(j+1))+1]*X[j+1]+
             1.*gField[4*((i-1)*Width+(j+1))+1]*X[j+1];
        
        dy *= (Y[i]-Y[i-1])*(X[j+1]-X[j])/64.0;
        dx *= (Y[i]-Y[i-1])*(X[j+1]-X[j])/64.0;
        
        a0[k] += dx;
        b0[k] += dy;

        dA = 9.*((label[    i*Width  +  j]==k)?1.0:0.0)+
             3.*((label[(i-1)*Width  +  j]==k)?1.0:0.0)+
             3.*((label[    i*Width+(j+1)]==k)?1.0:0.0)+
             1.*((label[(i-1)*Width+(j+1)]==k)?1.0:0.0);

        A[k] += dA*(Y[i]-Y[i-1])*(X[j+1]-X[j])/64.0;
        
        /*************************************************/

        // +- quadrant

        dgradU[0][1] = 9.*gField[4*(    i*Width     +j  )+1]+
                       3.*gField[4*((i+1)*Width     +j  )+1]+
                       3.*gField[4*(    i*Width   +(j-1))+1]+
                       1.*gField[4*((i+1)*Width   +(j-1))+1];

        dgradU[1][0] = 9.*gField[4*(    i*Width     +j  )+2]+
                       3.*gField[4*((i+1)*Width     +j  )+2]+
                       3.*gField[4*(    i*Width   +(j-1))+2]+
                       1.*gField[4*((i+1)*Width   +(j-1))+2];

        dgradU[0][1] *= (Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;
        dgradU[1][0] *= (Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;

        w[k] += dgradU[1][0]-dgradU[0][1];
        
        dy = 9.*gField[4*(    i*Width  +j  )+2]*Y[i]  +
             3.*gField[4*((i+1)*Width  +j  )+2]*Y[i+1]+
             3.*gField[4*(    i*Width+(j-1))+2]*Y[i]  +
             1.*gField[4*((i+1)*Width+(j-1))+2]*Y[i+1];

        dy-= 9.*gField[4*(    i*Width  +j  )+1]*Y[i]  +
             3.*gField[4*((i+1)*Width  +j  )+1]*Y[i+1]+
             3.*gField[4*(    i*Width+(j-1))+1]*Y[i]  +
             1.*gField[4*((i+1)*Width+(j-1))+1]*Y[i+1];
        
        dx = 9.*gField[4*(    i*Width  +j  )+2]*X[j]  +
             3.*gField[4*((i+1)*Width  +j  )+2]*X[j]  +
             3.*gField[4*(    i*Width+(j-1))+2]*X[j-1]+
             1.*gField[4*((i+1)*Width+(j-1))+2]*X[j-1];

        dx-= 9.*gField[4*(    i*Width  +j  )+1]*X[j]  +
             3.*gField[4*((i+1)*Width  +j  )+1]*X[j]  +
             3.*gField[4*(    i*Width+(j-1))+1]*X[j-1]+
             1.*gField[4*((i+1)*Width+(j-1))+1]*X[j-1];
        
        dy *= (Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;
        dx *= (Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;
        
        a0[k] += dx;
        b0[k] += dy;

        dA = 9.*((label[    i*Width  +  j]==k)?1.0:0.0)+
             3.*((label[(i+1)*Width  +  j]==k)?1.0:0.0)+
             3.*((label[    i*Width+(j-1)]==k)?1.0:0.0)+
             1.*((label[(i+1)*Width+(j-1)]==k)?1.0:0.0);

        A[k] += dA*(Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;
        
        /*************************************************/
        
        // -- quadrant

        dgradU[0][1] = 9.*gField[4*(    i*Width     +j  )+1]+
                       3.*gField[4*((i-1)*Width     +j  )+1]+
                       3.*gField[4*(    i*Width   +(j-1))+1]+
                       1.*gField[4*((i-1)*Width   +(j-1))+1];

        dgradU[1][0] = 9.*gField[4*(    i*Width     +j  )+2]+
                       3.*gField[4*((i-1)*Width     +j  )+2]+
                       3.*gField[4*(    i*Width   +(j-1))+2]+
                       1.*gField[4*((i-1)*Width   +(j-1))+2];

        dgradU[0][1] *= (Y[i]-Y[i-1])*(X[j]-X[j-1])/64.0;
        dgradU[1][0] *= (Y[i]-Y[i-1])*(X[j]-X[j-1])/64.0;

        w[k] += dgradU[1][0]-dgradU[0][1];
        
        dy = 9.*gField[4*(    i*Width  +j  )+2]*Y[i]  +
             3.*gField[4*((i-1)*Width  +j  )+2]*Y[i-1]+
             3.*gField[4*(    i*Width+(j-1))+2]*Y[i]  +
             1.*gField[4*((i-1)*Width+(j-1))+2]*Y[i-1];

        dy-= 9.*gField[4*(    i*Width  +j  )+1]*Y[i]  +
             3.*gField[4*((i-1)*Width  +j  )+1]*Y[i-1]+
             3.*gField[4*(    i*Width+(j-1))+1]*Y[i]  +
             1.*gField[4*((i-1)*Width+(j-1))+1]*Y[i-1];
        
        dx = 9.*gField[4*(    i*Width  +j  )+2]*X[j]  +
             3.*gField[4*((i-1)*Width  +j  )+2]*X[j]  +
             3.*gField[4*(    i*Width+(j-1))+2]*X[j-1]+
             1.*gField[4*((i-1)*Width+(j-1))+2]*X[j-1];

        dx-= 9.*gField[4*(    i*Width  +j  )+1]*X[j]  +
             3.*gField[4*((i-1)*Width  +j  )+1]*X[j]  +
             3.*gField[4*(    i*Width+(j-1))+1]*X[j-1]+
             1.*gField[4*((i-1)*Width+(j-1))+1]*X[j-1];
        
        dy *= (Y[i]-Y[i-1])*(X[j]-X[j-1])/64.0;
        dx *= (Y[i]-Y[i-1])*(X[j]-X[j-1])/64.0;
        
        a0[k] += dx;
        b0[k] += dy;

        dA = 9.*((label[    i*Width  +  j]==k)?1.0:0.0)+
             3.*((label[(i-1)*Width  +  j]==k)?1.0:0.0)+
             3.*((label[    i*Width+(j-1)]==k)?1.0:0.0)+
             1.*((label[(i-1)*Width+(j-1)]==k)?1.0:0.0);

        A[k] += dA*(Y[i+1]-Y[i])*(X[j]-X[j-1])/64.0;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)*sqrtf(2.); // Constant comming from lamb-oseen vortex;
    if(w[k]>0.){
      a=a0[k]/w[k]; 
      b=b0[k]/w[k];
    }
    else{
      a=0.;
      b=0.;
    } 
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