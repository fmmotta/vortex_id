#include <math.h>
#include "stencilManip.h"

int checkFourStencil(int i,int j,int Width,int Height){
  int neighbours=0;

  // 4-way stencil check
  if(bound_check(i+1,j,Width,Height))
    neighbours+=1;
  
  if(bound_check(i,j+1,Width,Height))
    neighbours+=1;

  if(bound_check(i-1,j,Width,Height))
    neighbours+=1;
  
  if(bound_check(i,j-1,Width,Height))
    neighbours+=1;
  
  return neighbours;
}

int updateFrom2DVelocityField(int Width,int Height,
                              float dx,float dy,
                              float *uField,float *lambField){
  int neighbours,i,j;
  float gradU[4],detU,trU,lamb2;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      neighbours = checkFourStencil(i,j,Width,Height);
      if(neighbours!=4)
        continue;
      
      // using centered differences to calculate the gradient
      gradU[0*2+0] = (uField[ux(i+1,j,Width)]-uField[ux(i-1,j,Width)])/(2.*dx);
      gradU[0*2+1] = (uField[uy(i+1,j,Width)]-uField[uy(i-1,j,Width)])/(2.*dx);
      gradU[1*2+0] = (uField[ux(i,j+1,Width)]-uField[ux(i,j-1,Width)])/(2.*dy);
      gradU[1*2+1] = (uField[uy(i,j+1,Width)]-uField[uy(i,j+1,Width)])/(2.*dy);
    
      detU = gradU[0]*gradU[3]-gradU[1]*gradU[2];
      trU = gradU[0]+gradU[3];

      // \Delta = (trU)^2 - 4.*detU ; \Delta<0 ==> Imaginary eigenvalues
      lamb2 = detU-(trU*trU)/4.;
      if(lamb2>=0)
        lambField[i*Width+j] = sqrt(lamb2);
      else
        lambField[i*Width+j] = 0.; // I might have to change to -1
    }

  return 0;
}

int uFieldToGradU(int Height,int Width,float *dx,
                  float *uField, float *gField){
  const int stcW=3,stcR=1;   // Stencil Width and Radius
  const int ouTW=14,inTW=ouTW+stcW-1; // Output Tile and Input Tile Width 
  int i,j,ii,jj,idx,jdx;
  int iBlks,jBlks;
  float uTile[inTW][inTW][2],gTile[ouTW][ouTW][4];
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;

  iBlks = (int) ((Height-1)/ouTW + 1);
  jBlks = (int) ( (Width-1)/ouTW + 1);

  for(i=0;i<iBlks;i+=1)
    for(j=0;j<jBlks;j+=1){
      
      // Loading the Input Tile
      for(ii=0;ii<inTw;ii+=1)
        for(jj=0;jj<inTw;jj+=1){
          idx = i*inTw + ii - stcR;
          jdx = i*inTw + jj - stcR;

          if( (idx<0) || (idx >= Height) ||
              (jdx<0) || (jdx >= Width) ){ 
            uTile[idx][jdx][0] = 0.;
            uTile[idx][jdx][1] = 0.;
          }
          else{
            uTile[ii][jj][0] = uField[2*(idx*Width+jdx)+0];
            uTile[ii][jj][1] = uField[2*(idx*Width+jdx)+1];
          }
        }

      // Specific finite diferences calculation
      // Rewrite this for general Stencil Calulations
      // Check latter the bounds of this loops
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          gTile[ii][jj][0] = (uTile[ii+1+1][jj+1][0] - 
                              uTile[ii+1-1][jj+1][0])/(2.*dx[0]);

          gTile[ii][jj][1] = (uTile[ii+1][jj+1+1][0] - 
                              uTile[ii+1][jj+1-1][0])/(2.*dx[1]);

          gTile[ii][jj][2] = (uTile[ii+1+1][jj+1][1] - 
                              uTile[ii+1-1][jj+1][1])/(2.*dx[0]);

          gTile[ii][jj][3] = (uTile[ii+1][jj+1+1][1] - 
                              uTile[ii+1][jj+1+1][1])/(2.*dx[1]);
        }
      
      // Writing the Output Tile
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          idx = i*ouTw + ii;
          jdx = i*ouTw + jj;

          if( (idx<=Height) && (jdx<=Width) ){
            gField[4*(idx*Width+jdx)+0] = gTile[ii][jj][0];
            gField[4*(idx*Width+jdx)+1] = gTile[ii][jj][1];
            gField[4*(idx*Width+jdx)+2] = gTile[ii][jj][2];
            gField[4*(idx*Width+jdx)+3] = gTile[ii][jj][3];
          }
        }
    }

  return 0;  
}

int multiChStencil(int Height,int Width,int Cha1,int Cha2,float *dx,
                   float *uField, float *gField){
  const int stcW=3,stcR=1;   // Stencil Width and Radius
  const int ouTW=14,inTW=ouTW+stcW-1; // Output Tile and Input Tile Width 
  int i,j,ii,jj,idx,jdx;
  int iBlks,jBlks;
  float uTile[inTW][inTW][2],gTile[ouTW][ouTW][4];
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;

  iBlks = (int) ((Height-1)/ouTW + 1);
  jBlks = (int) ( (Width-1)/ouTW + 1);

  for(i=0;i<iBlks;i+=1)
    for(j=0;j<jBlks;j+=1){
      
      // Loading the Input Tile
      for(ii=0;ii<inTw;ii+=1)
        for(jj=0;jj<inTw;jj+=1){
          idx = i*inTw + ii - stcR;
          jdx = i*inTw + jj - stcR;

          if( (idx<0) || (idx >= Height) ||
              (jdx<0) || (jdx >= Width) ){ 
            uTile[idx][jdx][0] = 0.;
            uTile[idx][jdx][1] = 0.;
          }
          else{
            uTile[ii][jj][0] = uField[2*(idx*Width+jdx)+0];
            uTile[ii][jj][1] = uField[2*(idx*Width+jdx)+1];
          }
        }

      // Specific finite diferences calculation
      // Rewrite this for general Stencil Calulations
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          gTile[ii][jj][0] = (uTile[ii+1+1][jj+1][0] - 
                              uTile[ii+1-1][jj+1][0])/(2.*dx[0]);

          gTile[ii][jj][1] = (uTile[ii+1][jj+1+1][0] - 
                              uTile[ii+1][jj+1-1][0])/(2.*dx[1]);

          gTile[ii][jj][2] = (uTile[ii+1+1][jj+1][1] - 
                              uTile[ii+1-1][jj+1][1])/(2.*dx[0]);

          gTile[ii][jj][3] = (uTile[ii+1][jj+1+1][1] - 
                              uTile[ii+1][jj+1+1][1])/(2.*dx[1]);
        }
      
      // Writing the Output Tile
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          idx = i*ouTw + ii;
          jdx = i*ouTw + jj;

          if( (idx<=Height) && (jdx<=Width) ){
            gField[4*(idx*Width+jdx)+0] = gTile[ii][jj][0];
            gField[4*(idx*Width+jdx)+1] = gTile[ii][jj][1];
            gField[4*(idx*Width+jdx)+2] = gTile[ii][jj][2];
            gField[4*(idx*Width+jdx)+3] = gTile[ii][jj][3];
          }
        }
    }

  return 0;  
}

int applyMaskConvChnWise(int Height,int Width,int Cha,float *dx,
                         int stcW,float *Mask,float *uField){
  int inTW=16,ouTW=inTw-stcW+1; // Output Tile and Input Tile Width 
  int stcR = stcW/2;            // Stencil/Mask  Radius
  int i,j,ii,jj,l,m,idx,jdx,k;
  int iBlks,jBlks;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if( stcW%2 == 0)
    return -5;

  if(ouTW<0){
    inTW = 32;
    ouTW = ouTW=inTw-stcW+1;
    if(ouTW<0)
      return -6;
  }

  iBlks = (int) ((Height-1)/ouTW + 1);
  jBlks = (int) ( (Width-1)/ouTW + 1);

  for(i=0;i<iBlks;i+=1)
    for(j=0;j<jBlks;j+=1){
      float uTile[inTW][inTW][Cha];
      float oTile[ouTW][ouTW][Cha];

      // Loading the Input Tile
      for(ii=0;ii<inTw;ii+=1)
        for(jj=0;jj<inTw;jj+=1){
          idx = i*inTw + ii - stcR;
          jdx = i*inTw + jj - stcR;

          if( (idx<0) || (idx >= Height) ||
              (jdx<0) || (jdx >= Width) ){ 
            for(k=0;k<Cha;k+=1)
              uTile[idx][jdx][k] = 0.;
          }
          else{
            for(k=0;k<Cha;k+=1)
              uTile[ii][jj][k] = uField[Cha*(idx*Width+jdx)+k];
          }
        }

      // Specific finite diferences calculation
      // Rewrite this for general Stencil Calulations
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){

          for(k=0;k<Cha;k+=1)         
            oTile[ii][jj][k]=0.;
          
          // Verify this Loop
          for(l=0;l<stcW;l+=1)
            for(m=0;m<stcW;m+=1)
              for(k=0;k<Cha;k+=1)
                oTile[ii][jj][k] += Mask[l*stcW+m]*uTile[ii+l][jj+m][k];
        }
      
      // Writing the Output Tile
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          idx = i*ouTw + ii;
          jdx = i*ouTw + jj;

          if( (idx<=Height) && (jdx<=Width) ){
            
            for(k=0;k<Cha;k+=1)
              uTile[Cha*(idx*Width+jdx)+k] = oTile[ii][jj][k];
          }
        }
    }

  return 0;
}

int gradUToGradU2(int Height,int Width,float *dx,
                  float *gField, float *g2Field){
  const int stcW=3,stcR=1;   // Stencil Width and Radius
  const int ouTW=14,inTW=ouTW+stcW-1; // Output Tile and Input Tile Width 
  int i,j,ii,jj,idx,jdx;
  int iBlks,jBlks;
  float gTile[inTW][inTW][4],g2Tile[ouTW][ouTW][4];
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(gField==NULL)
    return -3;
  if(g2Field==NULL)
    return -4;

  iBlks = (int) ((Height-1)/ouTW + 1);
  jBlks = (int) ( (Width-1)/ouTW + 1);

  for(i=0;i<iBlks;i+=1)
    for(j=0;j<jBlks;j+=1){
      
      // Loading the Input Tile
      for(ii=0;ii<inTw;ii+=1)
        for(jj=0;jj<inTw;jj+=1){
          idx = i*inTw + ii - stcR;
          jdx = i*inTw + jj - stcR;

          if( (idx<0) || (idx >= Height) ||
              (jdx<0) || (jdx >= Width) ){ 
            gTile[idx][jdx][0] = 0.;
            gTile[idx][jdx][1] = 0.;
            gTile[idx][jdx][2] = 0.;
            gTile[idx][jdx][3] = 0.;
          }
          else{
            gTile[ii][jj][0] = gField[2*(idx*Width+jdx)+0];
            gTile[ii][jj][1] = gField[2*(idx*Width+jdx)+1];
            gTile[ii][jj][2] = gField[2*(idx*Width+jdx)+2];
            gTile[ii][jj][3] = gField[2*(idx*Width+jdx)+3];
          }
        }

      // Converting Gradient of the Velocity Field
      // To the Gradient of the associated velocity Field
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          
          g2Tile[ii][jj][1]=(gTile[ii+1][jj+1+1][2] + gTile[ii+1][jj+1-1][2]
                            -gTile[ii+1][jj+1+1][1] - gTile[ii+1][jj+1-1][2]
                           -2.*gTile[ii+1][jj+1][2] +2.*gTile[ii+1][jj+1][2]);
          g2Tile[ii][jj][1]=g2Tile[ii][jj][1]/(2.*dx[1]*dx[1]);

          g2Tile[ii][jj][2]=(gTile[ii+1+1][jj+1][2] + gTile[ii+1-1][jj+1][2]
                            -gTile[ii+1+1][jj+1][1] - gTile[ii+1-1][jj+1][2]
                            -2.*gTile[ii+1][jj+1][2] +2.*gTile[ii+1][jj+1][2]);
          g2Tile[ii][jj][2]=-g2Tile[ii][jj][1]/(2.*dx[0]*dx[0]);

          g2Tile[ii][jj][0]=(gTile[ii+1+1][jj+1][2] + gTile[ii+1][jj+1+1][2]
                            -gTile[ii+1-1][jj+1][1] - gTile[ii+1][jj+1+1][2]);
          g2Tile[ii][jj][0]=g2Tile[ii][jj][0]/(4.*dx[0]*dx[1]);

          g2Tile[ii][jj][3] = - g2Tile[ii][jj][0];
        }
      
      // Writing the Output Tile
      for(ii=0;ii<ouTw;ii+=1)
        for(jj=0;jj<ouTw;jj+=1){
          idx = i*ouTw + ii;
          jdx = i*ouTw + jj;

          if( (idx<=Height) && (jdx<=Width) ){
            gField[4*(idx*Width+jdx)+0] = gTile[ii][jj][0];
            gField[4*(idx*Width+jdx)+1] = gTile[ii][jj][1];
            gField[4*(idx*Width+jdx)+2] = gTile[ii][jj][2];
            gField[4*(idx*Width+jdx)+3] = gTile[ii][jj][3];
          }
        }
    }

  return 0;  
}