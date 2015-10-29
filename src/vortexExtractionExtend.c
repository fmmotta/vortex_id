#include "lambdaInit.h"
#include "floodFill.h"
#include "vortexExtraction.h"

int vortexExt2ndSwirl(int Height,int Width, int nCnect,
                      double *x0, double *dx,double *sField,
                      double *gField,int *label,double **vCatalogOut)
{
  int i,j,k;
  double G,a,b,x,y,rc,gradU[2][2]; // vorticity
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

        gradU[0][0] = gField[4*(i*Width+j)+0];
        gradU[0][1] = gField[4*(i*Width+j)+1];
        gradU[1][0] = gField[4*(i*Width+j)+2];
        gradU[1][1] = gField[4*(i*Width+j)+3];
        
        y = x0[0] + i*dx[0]; // change x and y
        x = x0[1] + j*dx[1]; // change x and y

        A[k] += dx[0]*dx[1];
        w[k] += ( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        a0[k] += x*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
        b0[k] += y*( gradU[1][0]-gradU[0][1] )*dx[0]*dx[1];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)*sqrtf(2.); // Constant comming from lamb-oseen vortex;
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