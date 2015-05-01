#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../floodFill.h"

  float sField0[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     0.,1.,1.,0.,0.,0.,0.,1.,0.,1.,
                     1.,1.,0.,1.,0.,0.,1.,1.,1.,1.,
                     0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,
                     1.,1.,1.,0.,0.,0.,0.,0.,1.,1.,
                     1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

  float sField1[] = {1.,0.,1.,0.,1.,0.,0.,0.,0.,0.,
                     1.,0.,1.,1.,0.,0.,0.,1.,0.,1.,
                     1.,0.,0.,1.,0.,0.,1.,1.,1.,1.,
                     0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
                     1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,
                     1.,0.,1.,0.,0.,1.,0.,0.,1.,1.,
                     1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     0.,0.,1.,0.,0.,0.,1.,0.,0.,0. };

  float sField2[] = {1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                     1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                     1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                     1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                     1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,1.,0.,1.,0.,0.,0. };
  
  float sField3[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                     1.,0.,1.,1.,0.,1.,0.,0.,0.,1.,
                     1.,0.,0.,1.,0.,0.,0.,0.,1.,1.,
                     0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                     0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
                     1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,
                     1.,0.,1.,0.,0.,1.,0.,0.,1.,1.,
                     1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     0.,0.,1.,0.,0.,0.,1.,0.,0.,0. };

  float sField4[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,1.,
                     1.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                     1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                     0.,1.,0.,1.,0.,1.,0.,0.,1.,1.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     1.,1.,1.,1.,1.,0.,0.,1.,1.,1.,
                     0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     0.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                     0.,1.,0.,1.,0.,1.,0.,0.,0.,1.,
                     0.,1.,0.,1.,0.,1.,0.,0.,0.,1. };

  float sField5[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,1.,
                     1.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                     1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                     0.,1.,0.,0.,0.,0.,1.,0.,1.,1.,
                     0.,0.,1.,0.,0.,1.,0.,0.,0.,0.,
                     1.,0.,1.,1.,1.,0.,0.,1.,1.,1.,
                     1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                     1.,0.,1.,0.,0.,0.,1.,0.,1.,1.,
                     0.,1.,0.,1.,0.,1.,0.,0.,0.,1.,
                     0.,1.,0.,1.,0.,1.,0.,0.,0.,1. };


int main(int argc,char **argv){
  const int Width = 10, Height = 10, Pop=10,nVortex=3;
  int i,j,err,ngbr,found;
  int nbList[8],label[Width*Height],eqList[Pop];
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float x,y,v0y0 = 0.05;
  
  err = floodFill(sField0,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField0[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }
  
  err = floodFill(sField1,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField1[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }
  
  err = floodFill(sField2,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField2[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }
  
  err = floodFill(sField3,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField3[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }
  
  err = floodFill(sField4,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField4[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }

  err = floodFill(sField5,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField5[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }

  return 0;
}