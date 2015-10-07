#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "stencilExtended.h"

int main(){  
  const int Width = 20, Height = 20, padWidth=2;
  int i,j,err;
  float *X,*Xbuff,*Y,*Ybuff;

  X=(float*)malloc((Width+2*padWidth)*sizeof(float));
  if(X==NULL)
  	printf("problem allocating X\n");

  Xbuff=(float*)malloc((Width+2*padWidth)*sizeof(float));
  if(Xbuff==NULL)
  	printf("problem allocating Xbuff\n");

  Y=(float*)malloc((Height+2*padWidth)*sizeof(float));
  if(Y==NULL)
  	printf("problem allocating Y\n");

  Ybuff=(float*)malloc((Height+2*padWidth)*sizeof(float));
  if(Ybuff==NULL)
  	printf("problem allocating Ybuff\n");
  
  for(j=0;j<Width;j+=1)
  	X[j] = j;

  for(i=0;i<Height;i+=1)
  	Y[i] = 1000+i;
  
  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
  	printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
  	printf("problem in XtoXbuff - Y\n");

  printf("__ __ ");
  for(i=0;i<Width;i+=1)
    printf("%2.0f ",X[i]);
  printf("__ __ \n");

  for(i=0;i<(Width+2*padWidth);i+=1)
  	printf("%2.0f ",Xbuff[i]);
  printf("\n");

  printf("____ ____ ");
  for(i=0;i<Height;i+=1)
    printf("%.0f ",Y[i]);
  printf("____ ____ \n");

  for(i=0;i<(Height+2*padWidth);i+=1)
    printf("%.0f ",Ybuff[i]);
  printf("\n");

  free(X);
  free(Xbuff);
  free(Y);
  free(Ybuff);

  return 0;
}