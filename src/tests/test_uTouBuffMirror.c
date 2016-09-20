#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fieldSmoothing.h"
/*
  int uFieldTouBuff(int Height,int Width,double *uField,double *uBuff,int padWidth)
 */

int print_u(int Height,int Width,double *uField,int padWidth,int comp){
  int i,j,k;

  for(k=0;k<padWidth;k+=1){
  	for(j=0;j<(Width+2*padWidth);j+=1)
  	  printf("____ ");
    printf("\n");
  }

  for(i=0;i<Height;i+=1){
  	for(k=0;k<padWidth;k+=1)
      printf("____ "); 

    for(j=0;j<Width;j+=1)
      printf("%4.0f ",uField[2*(i*Width+j)+comp]);

    for(k=0;k<padWidth;k+=1)
      printf("____ "); 

    printf("\n");
  }

  for(k=0;k<padWidth;k+=1){
  	for(j=0;j<(Width+2*padWidth);j+=1)
  	  printf("____ ");
    printf("\n");
  }

  return 0;
}

int main(){  
  const int Width = 20, Height = 20, padWidth=2;
  int i,j,err;
  double *uField,*uBuff;
  
  uField = (double*)malloc(2*Height*Width*sizeof(double));
  if(uField==NULL){
    printf("problem allocking memmory\n");
    return -1;
  }
  
  for(i=0;i<Height;i+=1)
  	for(j=0;j<Width;j+=1){
      uField[2*(i*Width+j)+0] = 1000+(i+1)*Width+(j+1);
      uField[2*(i*Width+j)+1] = 2000+(i+1)*Width+(j+1);
  	}

  uBuff = (double*)malloc(2*(Height+2*padWidth)*(Width+2*padWidth)*sizeof(double));
  if(uBuff==NULL){
  	printf("problem allocking memmory\n");
    return -1;
  }

  for(i=0;i<(Height+2*padWidth);i+=1)
  	for(j=0;j<(Width+2*padWidth);j+=1){
      uBuff[2*(i*Width+j)+0] = 0.;
      uBuff[2*(i*Width+j)+1] = 0.;
  	}
 
  err = uFieldTouBuffMirror(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
  	printf("problems in uFieldTouBuff\n");
  
  err=print_u(Height,Width,uField,padWidth,0);
  
  printf("\n");
  for(i=0;i<(Height+2*padWidth);i+=1){
  	for(j=0;j<(Width+2*padWidth);j+=1)
      printf("%4.0f ",uBuff[2*(i*(Width+2*padWidth)+j)+0]);
    printf("\n");
  }
  
  free(uField);
  free(uBuff);
  
  return 0;
}