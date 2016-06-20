#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

void vortexAdaptiveQuickSort(double *v,int nCnect,int size,
                             int (*cmp)(const double*,const double*)){
  int i,j,k;
  double p[size],t[size];
  if(nCnect<2)
    return;

  for(k=0;k<size;k+=1)
    p[k] = v[size*(nCnect/2)+k];

  for(i=0,j=nCnect-1;;i++,j--){
    while( cmp(v+size*i,p) )
      i++;
    while( cmp(p,v+size*j) )
      j--;
    if(i >= j)
      break;
    
    for(k=0;k<size;k+=1){
      t[k] = v[size*i+k];
      v[size*i+k] = v[size*j+k];
      v[size*j+k] = t[k];
    }
  }
  vortexAdaptiveQuickSort(v,i,size,cmp);
  vortexAdaptiveQuickSort(v+size*i,nCnect-i,size,cmp);
}

inline int add_dgradU(int Height,int Width, int i,int j,int ik,int jk,
                      double *gField,double dgradU[][2],double X[],double Y[])
{
  double da;

  dgradU[0][0]=dgradU[0][1]=dgradU[1][0]=dgradU[1][1]=0.;
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      dgradU[0][0] = 9.*gField[4*(     i*Width +  j  )+0];
      dgradU[0][1] = 9.*gField[4*(     i*Width +  j  )+1];
      dgradU[1][0] = 9.*gField[4*(     i*Width +  j  )+2];
      dgradU[1][1] = 9.*gField[4*(     i*Width +  j  )+3];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0;
    }
    else{
      dgradU[0][0] = 9.*gField[4*(     i*Width   +   j  )+0]+
                     3.*gField[4*(     i*Width   +(j+jk))+0];

      dgradU[0][1] = 9.*gField[4*(     i*Width   +   j  )+1]+
                     3.*gField[4*(     i*Width   +(j+jk))+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*(     i*Width +(j+jk))+2];

      dgradU[1][1] = 9.*gField[4*(     i*Width +   j  )+3]+
                     3.*gField[4*(     i*Width +(j+jk))+3];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0;
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      dgradU[0][0] = 9.*gField[4*(     i*Width     +  j )+0]+
                     3.*gField[4*((i+ik)*Width     +  j )+0];

      dgradU[0][1] = 9.*gField[4*(     i*Width     +  j )+1]+
                     3.*gField[4*((i+ik)*Width     +  j )+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*((i+ik)*Width +   j  )+2];

      dgradU[1][1] = 9.*gField[4*(     i*Width +   j  )+3]+
                     3.*gField[4*((i+ik)*Width +   j  )+3];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0;
    }
    else{
      dgradU[0][0] = 9.*gField[4*(     i*Width     +  j )+0]+
                     3.*gField[4*((i+ik)*Width     +  j )+0]+
                     3.*gField[4*(     i*Width   +(j+jk))+0]+
                     1.*gField[4*((i+ik)*Width   +(j+jk))+0];

      dgradU[0][1] = 9.*gField[4*(     i*Width     +  j )+1]+
                     3.*gField[4*((i+ik)*Width     +  j )+1]+
                     3.*gField[4*(     i*Width   +(j+jk))+1]+
                     1.*gField[4*((i+ik)*Width   +(j+jk))+1];

      dgradU[1][0] = 9.*gField[4*(     i*Width +   j  )+2]+
                     3.*gField[4*((i+ik)*Width +   j  )+2]+
                     3.*gField[4*(     i*Width +(j+jk))+2]+
                     1.*gField[4*((i+ik)*Width +(j+jk))+2];

      dgradU[1][1] = 9.*gField[4*(     i*Width +   j  )+3]+
                     3.*gField[4*((i+ik)*Width +   j  )+3]+
                     3.*gField[4*(     i*Width +(j+jk))+3]+
                     1.*gField[4*((i+ik)*Width +(j+jk))+3];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  dgradU[0][0] *= da;
  dgradU[0][1] *= da;
  dgradU[1][0] *= da;
  dgradU[1][1] *= da;

  return 0;
}

inline int add_dxdy(int Height,int Width, int i,int j,int ik,int jk,
                    double *gField,double *dx,double *dy,double X[],double Y[])
{
  double da;
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *dy = 9.*gField[4*(     i*Width  + j  )+2]*Y[i];

      *dy-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i];
        
      *dx = 9.*gField[4*(     i*Width  + j  )+2]*X[j];

      *dx-= 9.*gField[4*(     i*Width  + j  )+1]*X[j];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
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

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
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

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
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

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *dy *= da;
  *dx *= da;

  return 0;
}

inline int add_dA(int Height,int Width, int i,int j,int k,int ik,int jk,
                  double *gField,int *label,double *dA,double X[],double Y[])
{
  double da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0);

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[     i*Width+(j+jk)]==k)?1.0:0.0);

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[(i+ik)*Width  +  j ]==k)?1.0:0.0);

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *dA = 9.*((label[     i*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[(i+ik)*Width  +  j ]==k)?1.0:0.0)+
            3.*((label[     i*Width+(j+jk)]==k)?1.0:0.0)+
            1.*((label[(i+ik)*Width+(j+jk)]==k)?1.0:0.0);

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
    
       
  *dA *= da;

  return 0;
}

int vortexExtFromVortCurv(int Height,int Width, int nCnect,double *X,double *Y,
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
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
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

inline int add_dx2dy2(int Height,int Width, int i,int j,int ik,int jk,
                      double *gField,double *XX,double *XY,double *YX,
                      double *YY,double X[],double Y[])
{
  double da;
  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *XX = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*X[j   ];
      *XX-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*X[j   ];

      *XY = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*Y[i   ];
      *XY-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*Y[i   ];

      *YX = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*X[j   ];
      *YX-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*X[j   ];
        
      *YY = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*Y[i   ];
      *YY-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*Y[i   ];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *XX = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk]*X[j+jk];
      *XX-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk]*X[j+jk];

      *XY = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk]*Y[i   ];
      *XY-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk]*Y[i   ];

      *YX = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i   ]*X[j+jk];
      *YX-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i   ]*X[j+jk];

      *YY = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i   ]*Y[i   ];
      *YY-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i   ]*Y[i   ];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      *XX = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*X[j   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j   ]*X[j   ];
      *XX-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*X[j   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j   ]*X[j   ];

      *XY = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j   ]*Y[i+ik];
      *XY-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j   ]*Y[i+ik];

      *YX = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*X[j   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik]*X[j   ];
      *YX-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*X[j   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik]*X[j   ];

      *YY = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik]*Y[i+ik];
      *YY-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik]*Y[i+ik];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *XX = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk]*X[j+jk]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j   ]*X[j   ]+
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*X[j+jk]*X[j+jk];
      *XX-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk]*X[j+jk]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j   ]*X[j   ]+
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*X[j+jk]*X[j+jk];

      *XY = 9.*gField[4*(     i*Width  + j  )+2]*X[j   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*X[j+jk]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*X[j   ]*Y[i+ik]+
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*X[j+jk]*Y[i+ik];
      *XY-= 9.*gField[4*(     i*Width  + j  )+1]*X[j   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*X[j+jk]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*X[j   ]*Y[i+ik]+
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*X[j+jk]*Y[i+ik];

      *YX = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*X[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i   ]*X[i+jk]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik]*X[i   ]+
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*Y[i+ik]*X[i+jk];
      *YX-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*X[j   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i   ]*X[j+jk]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik]*X[j   ]+
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*Y[i+ik]*X[j+jk];

      *YY = 9.*gField[4*(     i*Width  + j  )+2]*Y[i   ]*X[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+2]*Y[i   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+2]*Y[i+ik]*Y[i+ik]+
            1.*gField[4*((i+ik)*Width+(j+jk))+2]*Y[i+ik]*Y[i+ik];
      *YY-= 9.*gField[4*(     i*Width  + j  )+1]*Y[i   ]*Y[i   ]+
            3.*gField[4*(     i*Width+(j+jk))+1]*Y[i   ]*Y[i   ]+
            3.*gField[4*((i+ik)*Width  + j  )+1]*Y[i+ik]*Y[i+ik]+
            1.*gField[4*((i+ik)*Width+(j+jk))+1]*Y[i+ik]*Y[i+ik];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *XX *= da;
  *XY *= da;
  *YX *= da;
  *YY *= da;

  return 0;
}

inline int add_dx2dy2w2(int Height,int Width, int i,int j,int ik,int jk,
                        double *gField,double *XX,double *XY,double *YX,
                        double *YY,double X[],double Y[])
{
  double w2=0.,da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *XX   = 9.*w2*X[j   ]*X[j   ];
      *XY   = 9.*w2*X[j   ]*Y[i   ];
      *YX   = 9.*w2*Y[i   ]*X[j   ];
      *YY   = 9.*w2*Y[i   ]*Y[i   ];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *XX  = 9.*w2*X[j   ]*X[j   ];
      *XY  = 9.*w2*X[j   ]*Y[i   ];
      *YX  = 9.*w2*Y[i   ]*X[j   ];
      *YY  = 9.*w2*Y[i   ]*Y[i   ];
      
      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *XX += 3.*w2*X[j+jk]*X[j+jk];
      *XY += 3.*w2*X[j+jk]*Y[i   ];
      *YX += 3.*w2*Y[i   ]*X[j+jk];
      *YY += 3.*w2*Y[i   ]*Y[i   ];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *XX   = 9.*w2*X[j   ]*X[j   ];
      *XY   = 9.*w2*X[j   ]*Y[i   ];
      *YX   = 9.*w2*Y[i   ]*X[j   ];
      *YY   = 9.*w2*Y[i   ]*Y[i   ];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *XX += 3.*w2*X[j   ]*X[j   ];
      *XY += 3.*w2*X[j   ]*Y[i+ik];
      *YX += 3.*w2*Y[i+ik]*X[j   ];
      *YY += 3.*w2*Y[i+ik]*Y[i+ik];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *XX   = 9.*w2*X[j   ]*X[j   ];
      *XY   = 9.*w2*X[j   ]*Y[i   ];
      *YX   = 9.*w2*Y[i   ]*X[j   ];
      *YY   = 9.*w2*Y[i   ]*Y[i   ];

      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *XX += 3.*w2*X[j+jk]*X[j+jk];
      *XY += 3.*w2*X[j+jk]*Y[i   ];
      *YX += 3.*w2*Y[i   ]*X[j+jk];
      *YY += 3.*w2*Y[i   ]*Y[i   ];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *XX += 3.*w2*X[j   ]*X[j   ];
      *XY += 3.*w2*X[j   ]*Y[i+ik];
      *YX += 3.*w2*Y[i+ik]*X[j   ];
      *YY += 3.*w2*Y[i+ik]*Y[i+ik];

      w2 = gField[4*((i+ik)*Width+(j+jk))+2] 
         - gField[4*((i+ik)*Width+(j+jk))+1];
      w2 = w2*w2;

      *XX += 1.*w2*X[j+jk]*X[j+jk];
      *XY += 1.*w2*X[j+jk]*Y[i+ik];
      *YX += 1.*w2*Y[i+ik]*X[j+jk];
      *YY += 1.*w2*Y[i+ik]*Y[i+ik];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *XX *= da;
  *XY *= da;
  *YX *= da;
  *YY *= da;

  return 0;
}

int extractSecondMoment(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,double *vCatalog,
                        double *vortSndMomMatrix)
{
  int i,j,k,err;
  double dgradU[2][2]; // vorticity
  double w[nCnect],SndMom[4*nCnect],XX,XY,YX,YY;
  
  if((Height<=0)||(Width<=0))
    return -1;
  
  for(k=0;k<nCnect;k+=1){
    w[k]=0.;SndMom[4*k+0]=0;SndMom[4*k+1]=0;
    SndMom[4*k+2]=0;SndMom[4*k+3]=0;
  }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant        
        err=add_dgradU(Height,Width,i,j,1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        //w[k] += dgradU[1][0]-dgradU[0][1];
        w[k] += (dgradU[1][0]-dgradU[0][1])*(dgradU[1][0]-dgradU[0][1]);
        
        //XX = XY = YX = YY = 0.;
        //err=add_dx2dy2(Height,Width,i,j,1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        err=add_dx2dy2w2(Height,Width,i,j,1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        //w[k] += dgradU[1][0]-dgradU[0][1];
        w[k] += (dgradU[1][0]-dgradU[0][1])*(dgradU[1][0]-dgradU[0][1]);

        //XX = XY = YX = YY = 0.;
        //err=add_dx2dy2(Height,Width,i,j,-1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        err=add_dx2dy2w2(Height,Width,i,j,-1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        //w[k] += dgradU[1][0]-dgradU[0][1];
        w[k] += (dgradU[1][0]-dgradU[0][1])*(dgradU[1][0]-dgradU[0][1]);
        
        //XX = XY = YX = YY = 0.;
        //err=add_dx2dy2(Height,Width,i,j,1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        err=add_dx2dy2w2(Height,Width,i,j,1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        //w[k] += dgradU[1][0]-dgradU[0][1];
        w[k] += (dgradU[1][0]-dgradU[0][1])*(dgradU[1][0]-dgradU[0][1]);
        
        //XX = XY = YX = YY = 0.;
        //err=add_dx2dy2(Height,Width,i,j,-1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        err=add_dx2dy2w2(Height,Width,i,j,-1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    
    if(fabs(w[k])>0.){
      XX = SndMom[4*k+0]/w[k];
      XY = SndMom[4*k+1]/w[k];
      YX = SndMom[4*k+2]/w[k];
      YY = SndMom[4*k+3]/w[k];
    }
    else{
      XX = -1.;
      XY =  0.;
      YX =  0.;
      YY = -1.;
    }
    
    vortSndMomMatrix[4*k+0] = XX;
    vortSndMomMatrix[4*k+1] = XY;
    vortSndMomMatrix[4*k+2] = YX;
    vortSndMomMatrix[4*k+3] = YY;
  }

  return 0;
}

inline int add_dxdyw2(int Height,int Width, int i,int j,int ik,int jk,
                      double *gField,double *dx,double *dy,double X[],double Y[])
{
  double w2,da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *dy  = 9.*w2*Y[i   ];
      *dx  = 9.*w2*X[j   ];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *dy  = 9.*w2*Y[i   ];
      *dx  = 9.*w2*X[j   ];
      
      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *dy += 3.*w2*Y[i   ];
      *dx += 3.*w2*X[j+jk];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *dy  = 9.*w2*Y[i   ];
      *dx  = 9.*w2*X[j   ];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *dy += 3.*w2*Y[i+ik];
      *dx += 3.*w2*X[j   ];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *dy  = 9.*w2*Y[i   ];
      *dx  = 9.*w2*X[j   ];

      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *dy += 3.*w2*Y[i   ];
      *dx += 3.*w2*X[j+jk];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *dy += 3.*w2*Y[i+ik];
      *dx += 3.*w2*X[j   ];

      w2 = gField[4*((i+ik)*Width+(j+jk))+2] 
         - gField[4*((i+ik)*Width+(j+jk))+1];
      w2 = w2*w2;

      *dy += 1.*w2*Y[i+ik];
      *dx += 1.*w2*X[j+jk];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *dy *= da;
  *dx *= da;

  return 0;
}

inline int add_w2(int Height,int Width, int i,int j,int ik,int jk,
                  double *gField,double *w2out,double X[],double Y[])
{
  double w2,da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *w2out  = 9.*w2;

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *w2out  = 9.*w2;
      
      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *w2out += 3.*w2;

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *w2out  = 9.*w2;
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *w2out += 3.*w2;

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *w2out = 9.*w2;

      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *w2out += 3.*w2;
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *w2out += 3.*w2;

      w2 = gField[4*((i+ik)*Width+(j+jk))+2] 
         - gField[4*((i+ik)*Width+(j+jk))+1];
      w2 = w2*w2;

      *w2out += 1.*w2;

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *w2out *= da;

  return 0;
}

inline int add_wBkg(int Height,int Width, int i,int j,int ik,int jk,
                    double *wBkg,double *wBkgOut,double X[],double Y[])
{
  double da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      *wBkgOut = 9*wBkg[i*Width+j];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *wBkgOut = 9*wBkg[     i*Width  + j  ]
               + 3*wBkg[     i*Width+(j+jk)];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      *wBkgOut = 9*wBkg[     i*Width  + j  ]
               + 3*wBkg[(i+ik)*Width  + j  ];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      *wBkgOut = 9*wBkg[     i*Width  + j  ]
               + 3*wBkg[     i*Width+(j+jk)]
               + 3*wBkg[(i+ik)*Width  + j  ]
               + 1*wBkg[(i+ik)*Width+(j+jk)];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *wBkgOut *= da;

  return 0;
}

int extract012Momentsw2(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,
                        double *vCatalog,double *vortSndMomMatrix,
                        double *avgGradU)
{
  int i,j,k,err;
  double rc,a,b,G,dx,dy,dA,XX,XY,YX,YY,dgradU[2][2],dw2;
  double A[nCnect],a0[nCnect],b0[nCnect],w[nCnect]; // vorticity
  double SndMom[4*nCnect],w2[nCnect],gradU[4*nCnect];
  
  if((Height<=0)||(Width<=0))
    return -1;
  
  dgradU[0][0]=dgradU[0][1]=dgradU[1][0]=dgradU[1][1]=0.;
  for(k=0;k<nCnect;k+=1){
    w[k]=0.;A[k]=0.;a0[k]=0.;b0[k]=0.;w2[k]=0.;
    SndMom[4*k+0]=0.;SndMom[4*k+1]=0.;SndMom[4*k+2]=0.;SndMom[4*k+3]=0.;
    gradU[4*k+0]=0.;gradU[4*k+1]=0.;gradU[4*k+2]=0.;gradU[4*k+3]=0.;
  }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant        
        err=add_dgradU(Height,Width,i,j,1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        gradU[4*k+0] += dgradU[0][0];
        gradU[4*k+1] += dgradU[0][1];
        gradU[4*k+2] += dgradU[1][0];
        gradU[4*k+3] += dgradU[1][1];
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;
        
        err=add_dx2dy2w2(Height,Width,i,j,1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 

        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        gradU[4*k+0] += dgradU[0][0];
        gradU[4*k+1] += dgradU[0][1];
        gradU[4*k+2] += dgradU[1][0];
        gradU[4*k+3] += dgradU[1][1];
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,-1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dx2dy2w2(Height,Width,i,j,-1,1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        gradU[4*k+0] += dgradU[0][0];
        gradU[4*k+1] += dgradU[0][1];
        gradU[4*k+2] += dgradU[1][0];
        gradU[4*k+3] += dgradU[1][1];
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;
        
        err=add_dx2dy2w2(Height,Width,i,j,1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 

        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        gradU[4*k+0] += dgradU[0][0];
        gradU[4*k+1] += dgradU[0][1];
        gradU[4*k+2] += dgradU[1][0];
        gradU[4*k+3] += dgradU[1][1];
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,-1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,-1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;
        
        err=add_dx2dy2w2(Height,Width,i,j,-1,-1,gField,&XX,&XY,&YX,&YY,X,Y);
        if(err!=0) return err;
        SndMom[4*k+0] += XX; SndMom[4*k+1] += XY;
        SndMom[4*k+2] += YX; SndMom[4*k+3] += YY; 
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    
    rc= sqrt(A[k]/M_PI);

    if(fabs(w2[k])>0.){
      XX = SndMom[4*k+0]/w2[k];
      XY = SndMom[4*k+1]/w2[k];
      YX = SndMom[4*k+2]/w2[k];
      YY = SndMom[4*k+3]/w2[k];
      a=a0[k]/w2[k]; 
      b=b0[k]/w2[k];
    }
    else{
      XX = -1.;
      XY =  0.;
      YX =  0.;
      YY = -1.;
      a=X[0];
      b=Y[0];
    }
    
    G = w[k];

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;

    vortSndMomMatrix[4*k+0] = XX;
    vortSndMomMatrix[4*k+1] = XY;
    vortSndMomMatrix[4*k+2] = YX;
    vortSndMomMatrix[4*k+3] = YY;

    avgGradU[4*k+0] = gradU[4*k+0];
    avgGradU[4*k+1] = gradU[4*k+1];
    avgGradU[4*k+2] = gradU[4*k+2];
    avgGradU[4*k+3] = gradU[4*k+3];
  }

  return 0;
}

int extractAvgBkgVort(int Height,int Width,double *X,double *Y,
                      int nCnect,int *label,double *wBkg,
                      double *wBkgGammaOut)
{
  int i,j,k,err;
  double wBkgGamma[nCnect],dwBkgGamma;
  
  if((Height<=0)||(Width<=0))
    return -1;

  for(k=0;k<nCnect;k+=1)
    wBkgGamma[k] = 0.;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant        

        err=add_wBkg(Height,Width,i,j,1,1,wBkg,&dwBkgGamma,X,Y);
        if(err!=0) return err;
        wBkgGamma[k] += dwBkgGamma;

        /*************************************************/
        
        // -+ quadrant

        err=add_wBkg(Height,Width,i,j,-1,1,wBkg,&dwBkgGamma,X,Y);
        if(err!=0) return err;
        wBkgGamma[k] += dwBkgGamma;
        
        /*************************************************/

        // +- quadrant

        err=add_wBkg(Height,Width,i,j,1,-1,wBkg,&dwBkgGamma,X,Y);
        if(err!=0) return err;
        wBkgGamma[k] += dwBkgGamma;
        
        /*************************************************/
        
        // -- quadrant

        err=add_wBkg(Height,Width,i,j,-1,-1,wBkg,&dwBkgGamma,X,Y);
        if(err!=0) return err;
        wBkgGamma[k] += dwBkgGamma;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1)
    wBkgGammaOut[k] = wBkgGamma[k];

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
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;
        
        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];
        
        err=add_dxdy(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
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

int extLambOseenParams(int Height,int Width, int nCnect,double *X,double *Y,
                       double *sField,double *gField,int *label,
                       double *vCatalog)
{
  int i,j,k,err;
  double rc,a,b,G,dx,dy,dA,dgradU[2][2],dw2;
  double A[nCnect],a0[nCnect],b0[nCnect],w[nCnect],w2[nCnect]; 
  
  if((Height<=0)||(Width<=0))
    return -1;
  
  for(k=0;k<nCnect;k+=1){ w[k]=0.;A[k]=0.;a0[k]=0.;b0[k]=0.;w2[k]=0.; }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant        
        err=add_dgradU(Height,Width,i,j,1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        /*************************************************/
        
        // -+ quadrant
        err=add_dgradU(Height,Width,i,j,-1,1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,-1,1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,-1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,-1,1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;
        
        /*************************************************/

        // +- quadrant
        err=add_dgradU(Height,Width,i,j,1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;

        /*************************************************/
        
        // -- quadrant
        err=add_dgradU(Height,Width,i,j,-1,-1,gField,dgradU,X,Y);
        if(err!=0) return err;
        w[k] += dgradU[1][0]-dgradU[0][1];

        err=add_dA(Height,Width,i,j,k,-1,-1,gField,label,&dA,X,Y);
        if(err!=0) return err;
        A[k] += dA;

        err=add_w2(Height,Width,i,j,-1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;
        
        err=add_dxdyw2(Height,Width,i,j,-1,-1,gField,&dx,&dy,X,Y);
        if(err!=0) return err;
        a0[k] += dx;
        b0[k] += dy;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){
    
    rc= sqrt(A[k]/M_PI);

    if(fabs(w2[k])>0.){
      a=a0[k]/w2[k]; 
      b=b0[k]/w2[k];
    }
    else{
      a=X[0];
      b=Y[0];
    }
    
    G = w[k];

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  return 0;
}

inline int add_UVw2(int Height,int Width, int i,int j,int ik,int jk,
                    double *uField,double *gField,double *u,double *v,
                    double X[],double Y[])
{
  double w2,da;

  if((i+ik<0) || (i+ik)>=Height){
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *u  = 9.*w2*uField[2*(i*Width+j)+0];
      *v  = 9.*w2*uField[2*(i*Width+j)+1];

      da = fabs( (Y[i-ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *u  = 9.*w2*uField[2*(     i*Width  + j  )+0];
      *v  = 9.*w2*uField[2*(     i*Width  + j  )+1];
      
      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *u += 3.*w2*uField[2*(     i*Width+(j+jk))+0];
      *v += 3.*w2*uField[2*(     i*Width+(j+jk))+1];

      da = fabs( (Y[i-ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }
  else{
    if( (j+jk<0) || (j+jk)>=Width ){
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *u  = 9.*w2*uField[2*(     i*Width  + j  )+0];
      *v  = 9.*w2*uField[2*(     i*Width  + j  )+1];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *u += 3.*w2*uField[2*((i+ik)*Width  + j  )+0];
      *v += 3.*w2*uField[2*((i+ik)*Width  + j  )+1];

      da = fabs( (Y[i+ik]-Y[i])*(X[j-jk]-X[j]) )/64.0; 
    }
    else{
      w2 = gField[4*(     i*Width  + j  )+2] 
         - gField[4*(     i*Width  + j  )+1];
      w2 = w2*w2;

      *u  = 9.*w2*uField[2*(     i*Width  + j  )+0];
      *v  = 9.*w2*uField[2*(     i*Width  + j  )+1];

      w2 = gField[4*(     i*Width+(j+jk))+2] 
         - gField[4*(     i*Width+(j+jk))+1];
      w2 = w2*w2;

      *u += 3.*w2*uField[2*(     i*Width+(j+jk))+0];
      *v += 3.*w2*uField[2*(     i*Width+(j+jk))+1];
      
      w2 = gField[4*((i+ik)*Width  + j  )+2] 
         - gField[4*((i+ik)*Width  + j  )+1];
      w2 = w2*w2;

      *u += 3.*w2*uField[2*((i+ik)*Width  + j  )+0];
      *v += 3.*w2*uField[2*((i+ik)*Width  + j  )+1];

      w2 = gField[4*((i+ik)*Width+(j+jk))+2] 
         - gField[4*((i+ik)*Width+(j+jk))+1];
      w2 = w2*w2;

      *u += 1.*w2*uField[2*((i+ik)*Width+(j+jk))+0];
      *v += 1.*w2*uField[2*((i+ik)*Width+(j+jk))+1];

      da = fabs( (Y[i+ik]-Y[i])*(X[j+jk]-X[j]) )/64.0; 
    }
  }

  *u *= da;
  *v *= da;

  return 0;
}

int extVortexVelocity(int Height,int Width, int nCnect,double *X,double *Y,
                      double *uField,double *sField,double *gField,int *label,
                      double *uVort)
{
  int i,j,k,err;
  double dw2,U,V;
  double w2[nCnect],u[2*nCnect];
  
  if((Height<=0)||(Width<=0))
    return -1;
  
  for(k=0;k<nCnect;k+=1){ w2[k]=0.; u[2*k+0]=0.; u[2*k+1]=0.; }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      if((k>=0)&&(k<nCnect)){
        // ++ quadrant        
        err=add_w2(Height,Width,i,j,1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k] += dw2;

        err=add_UVw2(Height,Width,i,j,1,1,uField,gField,&U,&V,X,Y);
        if(err!=0) return err;
        u[2*k+0] += U;
        u[2*k+1] += V;

        /*************************************************/
        
        // -+ quadrant
        err=add_w2(Height,Width,i,j,-1,1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;

        err=add_UVw2(Height,Width,i,j,-1,1,uField,gField,&U,&V,X,Y);
        if(err!=0) return err;
        u[2*k+0] += U;
        u[2*k+1] += V;
        
        /*************************************************/

        // +- quadrant
        err=add_w2(Height,Width,i,j,1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;

        err=add_UVw2(Height,Width,i,j,1,-1,uField,gField,&U,&V,X,Y);
        if(err!=0) return err;
        u[2*k+0] += U;
        u[2*k+1] += V;

        /*************************************************/
        
        // -- quadrant
        err=add_w2(Height,Width,i,j,-1,-1,gField,&(dw2),X,Y);
        if(err!=0) return err;
        w2[k]+=dw2;

        err=add_UVw2(Height,Width,i,j,-1,-1,uField,gField,&U,&V,X,Y);
        if(err!=0) return err;
        u[2*k+0] += U;
        u[2*k+1] += V;
        
        /*************************************************/
      }
    }

  for(k=0;k<nCnect;k+=1){

    if(fabs(w2[k])>0.){
      U=u[2*k+0]/w2[k]; 
      V=u[2*k+1]/w2[k];
    }
    else{
      U=0.0;
      V=0.0;
    }

    uVort[2*k+0] = U;
    uVort[2*k+1] = V;
  }

  return 0;
}