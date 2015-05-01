#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

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

#define ux(i,j,w) (((i)*(w)+(j))*2+0)
#define uy(i,j,w) (((i)*(w)+(j))*2+1)

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