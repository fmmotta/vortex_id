#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lambdaInit.h"

int initZero(int Height,int Width, float **gFieldOut){
  int i,j,k;
  float *gField;

  if((*gFieldOut) == NULL){
    gField = (float*)malloc(4*Height*Width*sizeof(float));
    if(gField == NULL)
      return 1;
  }
  else
    gField = *gFieldOut;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gField[4*(i*Width+j)+0] = 0.;
      gField[4*(i*Width+j)+1] = 0.;
      gField[4*(i*Width+j)+2] = 0.;
      gField[4*(i*Width+j)+3] = 0.;
    }

  *gFieldOut = gField;
  return 0;
}

int gradUtoLamb(int Height,int Width, float *gField,float **sFieldOut){
  int i,j,k;
  float gradU[2][2],*sField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  
  if((*sFieldOut) == NULL){
    sField = (float*)malloc(Height*Width*sizeof(float));
    if(sField == NULL)
      return 1;
  }
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gField[4*(i*Width+j)+0];
      gradU[0][1] = gField[4*(i*Width+j)+1];
      gradU[1][0] = gField[4*(i*Width+j)+2];
      gradU[1][1] = gField[4*(i*Width+j)+3];

      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (gradU[0][0]*gradU[1][1]-gradU[0][1]*gradU[1][0]);
      lamb-= ((gradU[0][0]+gradU[1][1])*(gradU[0][0]+gradU[1][1]))/4.;
      
      if(lamb>0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }
    
  *sFieldOut = sField;

  return 0;
}

int addSingleOseen(int nVortex,float *parVortex, float *x0, float *dx, 
                   int Height,int Width, float **gFieldOut){
  int i,j,k;
  float gradU[2][2],*gField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;

  if(*gFieldOut==NULL)
    return 1;
  gField = *gFieldOut;
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{          
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      gField[4*(i*Width+j)+0]+=gradU[0][0];
      gField[4*(i*Width+j)+1]+=gradU[0][1];
      gField[4*(i*Width+j)+2]+=gradU[1][0];
      gField[4*(i*Width+j)+3]+=gradU[1][1];
    }

  return 0;
}

int addConstXYShear(float *x0, float *dx,int Height,
                    int Width, float v0y0,float **gFieldOut){
  int i,j,k;
  float gradU[2][2],*gField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;

  if(*gFieldOut==NULL)
    return 1;
  gField = *gFieldOut;
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];

      gradU[0][1] += v0y0; // got to improove this later
      
      gField[4*(i*Width+j)+0]+=gradU[0][0];
      gField[4*(i*Width+j)+1]+=gradU[0][1];
      gField[4*(i*Width+j)+2]+=gradU[1][0];
      gField[4*(i*Width+j)+3]+=gradU[1][1];
    }

  return 0;
}

/*
 * Add gradU calculation to the init functions because
 * it is necessary to further calculations for vortex extraction
 */

int initLambOseen2D(int nVortex,float *parVortex,
                    float *x0, float *dx, int Height,int Width,
                    float **sFieldOut){
  int i,j,k;
  float gradU[2][2],*sField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb;
  float cutoff=0.001; // cutoff should be adjustable
  
  if((*sFieldOut) == NULL){
    sField = (float*)malloc(Height*Width*sizeof(float));
    if(sField == NULL)
      return 1;
  }
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{          
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (gradU[0][0]*gradU[1][1]-gradU[0][1]*gradU[1][0]);
      lamb-= ((gradU[0][0]+gradU[1][1])*(gradU[0][0]+gradU[1][1]))/4.;
      
      if(lamb>0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }

  *sFieldOut = sField;
  return 0;
}

int initOseenShear2D(int nVortex,float *parVortex,
                     float *x0, float *dx, int Height,int Width,
                     float v0y0, float **sFieldOut){
  int i,j,k;
  float gradU[2][2],*sField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  
  if((*sFieldOut) == NULL){
    sField = (float*)malloc(Height*Width*sizeof(float));
    if(sField == NULL)
      return 1;
  }
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      gradU[0][1] += v0y0; // got to improove this later
      
      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (gradU[0][0]*gradU[1][1]-gradU[0][1]*gradU[1][0]);
      lamb-= ((gradU[0][0]+gradU[1][1])*(gradU[0][0]+gradU[1][1]))/4.;
      
      if(lamb>0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }
    
  *sFieldOut = sField;
  return 0;
}
