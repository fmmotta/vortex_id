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

int uFieldToGradUopenFOAM(int Height,int Width,int type,float *x0,float *dx,
                          float *X,float *Y, float *uField, float *gField)
{
  const int refHeight;
  const int stcW=3,stcR=1;   // Stencil Width and Radius
  const int ouTW=14,inTW=ouTW+stcW-1; // Output Tile and Input Tile Width 
  int i,j,ii,jj,idx,jdx,ip,jp;
  int iBlks,jBlks;
  float uTile[inTW][inTW][2],gTile[ouTW][ouTW][4];
  float um,u0,up,vm,v0,vp,dYm[Height],dYp[Height];
  float cmin[Height],c0[Height],cmax[Height];
  float d0,dm,dmm,b1,b2;

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

  iBlks = (int) ((Height-1)/ouTW + 1);
  jBlks = (int) ( (Width-1)/ouTW + 1);
  
  //dYm[0] = (Y[0]-0.)/2.; // Hardcode: 3 point stencil with wall
  dYm[0] = (Y[1]-Y[0]);
  for(i=1;i<Height;i+=1)
    dYm[i] = Y[i]-Y[i-1];
  
  for(i=0;i<Height-1;i+=1)
    dYp[i] = (Y[i+1]-Y[i]);
  //dYp[Height-1] = 1.0-Y[Height-1];
  dYp[Height-1] = Y[Height-1]-Y[Height-2];

  printf("dx[0] = %f\n",dx[0]);
  // 1D, 3-point, 1st derivative stencil weights
  for(i=0;i<Height;i+=1){
    cmin[i] = -dYp[i]/(dYm[i]*(dYm[i]+dYp[i]));
      c0[i] =  (dYp[i]-dYm[i])/(dYm[i]*dYp[i]);
    cmax[i] =  dYm[i]/(dYp[i]*(dYm[i]+dYp[i]));
  }

  // 1D, separate stencil for the mid of the channel
  b1  = Y[Height-1] - Y[Height-2];
  b2  = Y[Height-1] - Y[Height-3];
  d0  =  (b1+b2)/(b1*b2);
  dm  = -b2/(b1*(b2-b1));
  dmm =  b1/(b2*(b2-b1));
  
  for(i=0;i<Height;i+=1){
    if(i==0){
      for(j=0;j<Width;j+=1){
        ip = i; jp = (j+Width-1)%Width;
        um = uField[ 2*(ip*Width+jp) + 0 ];
        vm = uField[ 2*(ip*Width+jp) + 1 ];
        
        u0 = uField[ 2*(i*Width+j) + 0 ];
        v0 = uField[ 2*(i*Width+j) + 1 ];
        
        ip = i; jp = (j+1)%Width;
        up = uField[ 2*(ip*Width+jp) + 0 ];
        vp = uField[ 2*(ip*Width+jp) + 1 ];

        // x-derivatives
        gField[4*(i*Width+j)+0] = (up-um)/(2.*dx[0]);
        gField[4*(i*Width+j)+2] = (vp-vm)/(2.*dx[0]);

        // ------
 
        ip = i-1; jp = j; // Down the wall
        um = 0.;          // No slip condition
        vm = 0.;          // No slip condition
        
        // just above the wall
        u0 = uField[ 2*(i*Width+j) + 0 ];
        v0 = uField[ 2*(i*Width+j) + 1 ];
        
        ip = i+1; jp = j;
        up = uField[ 2*(ip*Width+jp) + 0 ];
        vp = uField[ 2*(ip*Width+jp) + 1 ];

        // y-derivatives
        gField[4*(i*Width+j)+1] = cmin[i]*um + c0[i]*u0 + cmax[i]*up;
        gField[4*(i*Width+j)+3] = cmin[i]*vm + c0[i]*v0 + cmax[i]*vp;
      }
    }
    else if(i==Height-1){
      for(j=0;j<Width;j+=1){
        ip = i; jp = (j+Width-1)%Width;
        um = uField[ 2*(ip*Width+jp) + 0 ];
        vm = uField[ 2*(ip*Width+jp) + 1 ];
        
        u0 = uField[ 2*(i*Width+j) + 0 ];
        v0 = uField[ 2*(i*Width+j) + 1 ];
        
        ip = i; jp = (j+1)%Width;
        up = uField[ 2*(ip*Width+jp) + 0 ];
        vp = uField[ 2*(ip*Width+jp) + 1 ];

        // x-derivatives
        gField[4*(i*Width+j)+0] = (up-um)/(2.*dx[0]);
        gField[4*(i*Width+j)+2] = (vp-vm)/(2.*dx[0]);

        ip = i; jp = j;
        u0 = uField[ 2*(ip*Width+jp) + 0 ];
        v0 = uField[ 2*(ip*Width+jp) + 1 ];

        ip = i-1; jp = j;
        um = uField[ 2*(ip*Width+jp) + 0 ];
        vm = uField[ 2*(ip*Width+jp) + 1 ];

        ip = i-2; jp = j;
        up = uField[ 2*(ip*Width+jp) + 0 ];
        vp = uField[ 2*(ip*Width+jp) + 1 ];

        // y-derivatives
        gField[4*(i*Width+j)+1] = d0*u0 + dm*um + dmm*up;
        gField[4*(i*Width+j)+3] = d0*v0 + dm*vm + dmm*vp;        
      }
    }
    else{
      for(j=0;j<Width;j+=1){
        ip = i; jp = (j+Width-1)%Width;
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
        
        ip = i-1; jp = j; 
        um = uField[ 2*(ip*Width+jp) + 0 ];
        vm = uField[ 2*(ip*Width+jp) + 1 ];
        
        ip = i; jp = j; 
        u0 = uField[ 2*(ip*Width+jp) + 0 ];
        v0 = uField[ 2*(ip*Width+jp) + 1 ];
        
        ip = i+1; jp = j;
        up = uField[ 2*(ip*Width+jp) + 0 ];
        vp = uField[ 2*(ip*Width+jp) + 1 ];

        // y-derivatives
        gField[4*(i*Width+j)+1] = cmin[i]*um + c0[i]*u0 + cmax[i]*up;
        gField[4*(i*Width+j)+3] = cmin[i]*vm + c0[i]*v0 + cmax[i]*vp;
      }
    }
  }

  return 0;
}

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
      gField[4*(i*Width+j)+0] = (up-um)/(2.*dx[0]);
      gField[4*(i*Width+j)+2] = (vp-vm)/(2.*dx[0]);
      

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
      gField[4*(i*Width+j)+1] = (up-um)/(2.*dx[1]);
      gField[4*(i*Width+j)+3] = (vp-vm)/(2.*dx[1]);
    }

  return 0;
}

int UToGradnonUnifFrame(int Height,int Width,int type,float *x0,float *dx,
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

  for(i=1;i<Height-1;i+=1)
    for(j=1;j<Width-1;j+=1){
      ay = Y[i]-Y[i-1];
      by = Y[i+1]-Y[i];
      ax = X[j]-X[j-1];
      bx = X[j+1]-X[j];

      ip = (i-1+Height)%Height; jp = j;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = (i+1)%Height; jp = j;
      up = uField[ 2*(ip*Width+jp) + 0 ];
      vp = uField[ 2*(ip*Width+jp) + 1 ];

      cm = - bx/(ax*(ax+bx));
      c0 =   (bx-ax)/(ax*bx);
      cp =   ax/(bx*(ax+bx));

      // x-derivatives
      gField[4*(i*Width+j)+0] = cm*um + c0*u0 + cp*up;
      gField[4*(i*Width+j)+2] = cm*vm + c0*v0 + cp*vp;

      ip = i; jp = (j-1+Width)%Width;
      um = uField[ 2*(ip*Width+jp) + 0 ];
      vm = uField[ 2*(ip*Width+jp) + 1 ];
      
      ip = i; jp = j;
      u0 = uField[ 2*(ip*Width+jp) + 0 ];
      v0 = uField[ 2*(ip*Width+jp) + 1 ];

      ip = i; jp = (j+1)%Width;
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