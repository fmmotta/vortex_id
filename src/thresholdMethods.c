#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

int componentMaxima(int Height, int Width, int *label, double *sField,
	                int Nmax,double *cMax)
{
  int i,j;

  if(Height<=0 || Width<=0)
    return -1;

  if(sField==NULL)
    return -2;

  for(i=0;i<Nmax;i+=1)
  	cMax[i] = -1.;

  for(i=0;i<Height;i+=1)
  	for(j=0;j<Width;j+=1){
  	  if(label[i*Width+j]>=0){
        if(cMax[label[i*Width+j]]<sField[i*Width+j])
          cMax[label[i*Width+j]] = sField[i*Width+j];
  	  }
  	}

  return 0;
}

int applyGlobalThreshold(int Height,int Width,double *sField,double theta)
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
      if(sField[i*Width+j]<=theta)
        sField[i*Width+j]=0.;
    }
  
  return 0;
}

int applyComponentThreshold(int Height,int Width,int *label,
	                        double *sField,double *theta)
{
  int i,j;

  if(Height<=0 || Width<=0)
    return -1;

  if(sField==NULL)
    return -2;

  if(theta==NULL)
    return -3;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      if(label[i*Width+j]>=0 && sField[i*Width+j]<=theta[label[i*Width+j]] )
          sField[i*Width+j]=0.;
    }
  
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
  
  // modify this part
  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)*(sqrtf(2.)); // Constant comming from 
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
    rc= sqrt(A[k]/M_PI)/1.12091; // Constant comming from 
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
