#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int vortexExtraction(int Height,int Width, int nCnect,
	                 float *x0, float *dx,float *sField,
	                 float *gField,int *label,float **vCatalogOut)
{
  int i,j,k;
  float G,a,b,x,y,rc,gradU[2][2]; // vorticity
  float w[nCnect],A[nCnect],a0[nCnect],b0[nCnect];
  float vCatalog[4*nCnect];

  printf("nCnect=%d\n",nCnect);

  if((Height<=0)||(Width<=0))
  	return -1;

  if(*vCatalogOut!=NULL)
    return -2;

  for(i=0;i<nCnect;i+=1){
    w[k]=0.;A[k]=0.;a0[k]=0.;b0[k]=0.;
  }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      k=label[i*Width+j];

      gradU[0][0] = gField[4*(i*Width+j)+0];
      gradU[0][1] = gField[4*(i*Width+j)+1];
      gradU[1][0] = gField[4*(i*Width+j)+2];
      gradU[1][1] = gField[4*(i*Width+j)+3];

      if((k>=0)&&(k<nCnect)){
        x = x0[0] + j*dx[0];
        y = x0[1] + i*dx[1];

        A[k] += dx[0]*dx[1];
        //w[k] += fabs(gradU[0][1]-gradU[1][0])*dx[0]*dy[0];// maybe remove fabs
        w[k] += ((gradU[1][0]-gradU[0][1])/2.)*dx[0]*dx[1];
        //a0[k] += a*dx[0]*dx[1];
        //b0[k] += b*dx[0]*dx[1];
        a0[k] += x*w[k];
        b0[k] += y*w[k];
      }
    }

  for(k=0;k<nCnect;k+=1){
    rc= sqrt(A[k]/M_PI)/1.12091; // Constant comming from lamb-oseen vortex;
    a=a0[k]/w[k]; 
    b=b0[k]/w[k]; 
    //a = a0[k]/A[k];
    //b = b0[k]/A[k];
    G = (2.967284274)*w[k]; // Constant comming from lamb-oseen vortex;

    vCatalog[4*k+0] = G;
    vCatalog[4*k+1] = rc;
    vCatalog[4*k+2] = a;
    vCatalog[4*k+3] = b;
  }

  *vCatalogOut=vCatalog;

  return 0;
}