#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "stencilManip.h"

// Y is between 0 and 1 (half channel)
// It's allways periodic on X (can be stream or span-wise)
//  
//     |
//     |
//     |
//     |
//  /\ |
//  || |
//  || | 
//  i  |
//     |
//     |_ _ _ _ _ _ _ _ _ _ _ _ _ _ 
//       
//        j ==> 
//

int UToGradUuniformTorus(int Height,int Width,int type,float *x0,float *dx,
                         float *X,float *Y, float *uField, float *gField){
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){      

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      // x-derivatives
      gField[4*(i*Width+j)+0] = (up-um)/(2.*dx[0]);
      gField[4*(i*Width+j)+2] = (vp-vm)/(2.*dx[0]);

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      // x-derivatives
      gField[4*(i*Width+j)+1] = (up-um)/(2.*dx[1]);
      gField[4*(i*Width+j)+3] = (vp-vm)/(2.*dx[1]);
    }

  return 0;
}

int UToGradUnonUnifFrame(int Height,int Width,int type,float *x0,float *dx,
                         float *X,float *Y, float *uField, float *gField){
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;

  // Why this disposition works? it shouldn't work....
  // maybe is of that 'future' clauses that should be modified
  // result: changing the future clauses (x->y y->x) solved the 
  // problem and will be the default now. This comment should be
  // removed soon

  for(i=1;i<Height-1;i+=1)
    for(j=1;j<Width-1;j+=1){
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];
      
      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;
      
      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }

  return 0;
} 

int UToGradUnUnifHalfFrame(int Height,int Width,int type,float *x0,float *dx,
                           float *X,float *Y, float *uField, float *gField){
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;

  for(i=1;i<Height-1;i+=1){
    j=0;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      
      ax = X[j+1]-X[j];
      bx = X[j+2]-X[j];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = -(ax+bx)/(ax*bx);
      cm =  bx/(ax*(bx-ax));
      cp = -ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }

    for(j=1;j<Width-1;j+=1){
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }

    j=Width-1;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      
      ax = X[j]-X[j-1];
      bx = X[j]-X[j-2];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = (ax+bx)/(ax*bx);
      cm = - bx/(ax*(bx-ax));
      cp =  ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }
  }

  return 0;
}

int UToGradUnUnifFullFrame(int Height,int Width,int type,float *x0,float *dx,
                           float *X,float *Y, float *uField, float *gField){
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;
  
  i=0;
  {
    j=0;
    {
      ay = Y[i+1]-Y[i];
      by = Y[i+2]-Y[i];      
      ax = X[j+1]-X[j];
      bx = X[j+2]-X[j];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = -(ax+bx)/(ax*bx);
      cm =  bx/(ax*(bx-ax));
      cp = -ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 = -(ay+by)/(ay*by);
      cm =  by/(ay*(by-ay));
      cp = -ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }

    for(j=1;j<Width-1;j+=1){
      ay = Y[i+1]-Y[i];
      by = Y[i+2]-Y[i];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;

      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 = -(ay+by)/(ay*by);
      cm =  by/(ay*(by-ay));
      cp = -ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }

    j=Width-1;
    {
      ay = Y[i+1]-Y[i];
      by = Y[i+2]-Y[i];      
      ax = X[j]-X[j-1];
      bx = X[j]-X[j-2];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = (ax+bx)/(ax*bx);
      cm = - bx/(ax*(bx-ax));
      cp =  ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i+2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 = -(ay+by)/(ay*by);
      cm =  by/(ay*(by-ay));
      cp = -ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }
  }

  for(i=1;i<Height-1;i+=1){
    j=0;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      
      ax = X[j+1]-X[j];
      bx = X[j+2]-X[j];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = -(ax+bx)/(ax*bx);
      cm =  bx/(ax*(bx-ax));
      cp = -ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }

    for(j=1;j<Width-1;j+=1){
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }

    j=Width-1;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      
      ax = X[j]-X[j-1];
      bx = X[j]-X[j-2];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = (ax+bx)/(ax*bx);
      cm = - bx/(ax*(bx-ax));
      cp =  ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      cm = - by/(ay*(ay+by));
      c0 =   (by-ay)/(ay*by);
      cp =   ay/(by*(ay+by));

      // y-derivatives
      gField[4*(i*Width+j)+1] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+3] = cm*vm + c0*v0 + cp*vp;
    }
  }

  i=Height-1;
  {
    j=0;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i]-Y[i-2];      
      ax = X[j+1]-X[j];
      bx = X[j+2]-X[j];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j+2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = -(ax+bx)/(ax*bx);
      cm =  bx/(ax*(bx-ax));
      cp = -ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 =  (ay+by)/(ay*by);
      cm = -by/(ay*(by-ay));
      cp =  ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }

    for(j=1;j<Width-1;j+=1){
      ay = Y[i]-Y[i-1];
      by = Y[i]-Y[i-2];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 =  (ay+by)/(ay*by);
      cm = -by/(ay*(by-ay));
      cp =  ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }

    j=Width-1;
    {
      ay = Y[i]-Y[i-1];
      by = Y[i]-Y[i-2];      
      ax = X[j]-X[j-1];
      bx = X[j]-X[j-2];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-1;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = j-2;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      c0 = (ax+bx)/(ax*bx);
      cm = - bx/(ax*(bx-ax));
      cp =  ax/(bx*(bx-ax));

      // x-derivatives
      gField[4*(i*Width+j)+0] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+2] = c0*v0 + cm*vm + cp*vp;
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-1; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i-2; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];
      
      c0 =  (ay+by)/(ay*by);
      cm = -by/(ay*(by-ay));
      cp =  ay/(by*(by-ay));

      // y-derivatives
      gField[4*(i*Width+j)+1] = c0*u0 + cm*um + cp*up;
      gField[4*(i*Width+j)+3] = c0*v0 + cm*vm + cp*vp;
    }
  }

  return 0;
}