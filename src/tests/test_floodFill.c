#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"

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
  int nbList[8],label[Width*Height],eqList[Pop],**eqClass;
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float x,y,v0y0 = 0.05;

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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
  
  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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
  
  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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
  
  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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
  
  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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

  err = floodFill(sField,Width,Height,eqClass,label);
  err = renameLabels(Width,Height,label);

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

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  return 0;
}