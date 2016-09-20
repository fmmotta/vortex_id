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
  	  printf("______ ");
    printf("\n");
  }

  for(i=0;i<Height;i+=1){
  	for(k=0;k<padWidth;k+=1)
      printf("______ "); 

    for(j=0;j<Width;j+=1)
      printf("%2.3f ",uField[2*(i*Width+j)+comp]);

    for(k=0;k<padWidth;k+=1)
      printf("______ "); 

    printf("\n");
  }

  for(k=0;k<padWidth;k+=1){
  	for(j=0;j<(Width+2*padWidth);j+=1)
  	  printf("______ ");
    printf("\n");
  }

  return 0;
}

int main(){  
  const int Width = 20, Height = 20, padWidth=3;
  int i,j,err;
  double *X,*Xbuff,*Y,*Ybuff;
  double *uField,*uBuff,*uFilt;
  
  X=(double*)malloc((Width+2*padWidth)*sizeof(double));
  if(X==NULL)
    printf("problem allocating X\n");

  Xbuff=(double*)malloc((Width+2*padWidth)*sizeof(double));
  if(Xbuff==NULL)
    printf("problem allocating Xbuff\n");

  Y=(double*)malloc((Height+2*padWidth)*sizeof(double));
  if(Y==NULL)
    printf("problem allocating Y\n");

  Ybuff=(double*)malloc((Height+2*padWidth)*sizeof(double));
  if(Ybuff==NULL)
    printf("problem allocating Ybuff\n");

  for(j=0;j<Width;j+=1)
    X[j] = ((double)j)/10.;

  for(i=0;i<Height;i+=1)
    Y[i] = ((double)i)/10.;

  for(i=0;i<Height;i+=1)
    printf("%lf ",Y[i]);
  printf("\n");
  for(j=0;j<Width;j+=1)
    printf("%lf ",Y[j]);
  printf("\n");
  
  err = XtoXbuffMirror(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuffMirror(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  uField = (double*)malloc(2*Height*Width*sizeof(double));
  if(uField==NULL){
    printf("problem allocking memmory\n");
    return -1;
  }

  uFilt = (double*)malloc(2*Height*Width*sizeof(double));
  if(uFilt==NULL){
    printf("problem allocking memmory\n");
    return -1;
  }

  uBuff = (double*)malloc(2*(Height+2*padWidth)*(Width+2*padWidth)*sizeof(double));
  if(uBuff==NULL){
    printf("problem allocking memmory\n");
    return -1;
  }
  
  for(i=0;i<Height;i+=1)
  	for(j=0;j<Width;j+=1){
      if(j>(Width/2))
        uField[2*(i*Width+j)+0] = 1.;
      else
        uField[2*(i*Width+j)+0] = 0.;

      if(i>(Height/2))
        uField[2*(i*Width+j)+1] = -1.;
      else
        uField[2*(i*Width+j)+1] = 0.;
  	}

  for(i=0;i<(Height+2*padWidth);i+=1)
    for(j=0;j<(Width+2*padWidth);j+=1){
      uBuff[2*(i*Width+j)+0] = 0.;
      uBuff[2*(i*Width+j)+1] = 0.;
    }
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uFilt[2*(i*Width+j)+0] = 0.;
      uFilt[2*(i*Width+j)+1] = 0.;
    }
 
  err = uFieldTouBuffMirror(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
  	printf("problems in uFieldTouBuff\n");
  
  printf("un-buffed u:\n");
  err=print_u(Height,Width,uField,padWidth,0);
  printf("un-buffed v:\n");
  err=print_u(Height,Width,uField,padWidth,1);
  
  printf("Un-Filtered u:\n");
  printf("\n");
  for(i=0;i<(Height+2*padWidth);i+=1){
  	for(j=0;j<(Width+2*padWidth);j+=1)
      printf("%4.0f ",uBuff[2*(i*(Width+2*padWidth)+j)+0]);
    printf("\n");
  }

  printf("Un-Filtered v:\n");
  printf("\n");
  for(i=0;i<(Height+2*padWidth);i+=1){
    for(j=0;j<(Width+2*padWidth);j+=1)
      printf("%4.0f ",uBuff[2*(i*(Width+2*padWidth)+j)+1]);
    printf("\n");
  }

  //err=gaussianFilterNonUniform(Height,Width,padWidth,gauss5x5mask,Xbuff,Ybuff,uBuff,uFilt);

  //double sigma2 = (0.84089642*0.84089642)*(0.2*0.2)*(padWidth/2);
  double sigma2 = 0.84089642*0.84089642*0.1*0.1;
  err=gaussianFilterNonUniform(Height,Width,padWidth,Xbuff,Ybuff,uBuff,sigma2,uFilt);
  if(err!=0){
    printf("Problems in Gaussian non uniform\n");

    free(uField);
    free(uBuff);
    free(uFilt);
    free(X);
    free(Y);
    free(Xbuff);
    free(Ybuff);
    
    return -1;
  }

  printf("Filtered u:\n");
  printf("\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.3f ",uFilt[2*(i*Width+j)+0]);
    printf("\n");
  }

  printf("Filtered v:\n");
  printf("\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.3f ",uFilt[2*(i*Width+j)+1]);
    printf("\n");
  }

  printf("\n\n\n");

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uField[2*(i*Width+j)+0] = 0.;
      uField[2*(i*Width+j)+1] = 0.;
    }

  uField[2*((Height/2)*Width+(Width/2))+0] =  1.;
  uField[2*((Height/2)*Width+(Width/2))+1] = -1.;
  err = uFieldTouBuffMirror(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
    printf("problems in uFieldTouBuff\n");

  err=gaussianFilterNonUniform(Height,Width,padWidth,Xbuff,Ybuff,uBuff,sigma2,uFilt);
  if(err!=0){
    printf("Problems in Gaussian non uniform\n");

    free(uField);
    free(uBuff);
    free(uFilt);
    free(X);
    free(Y);
    free(Xbuff);
    free(Ybuff);
    
    return -1;
  }

  printf("Filtered u:\n");
  printf("\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.7f ",uFilt[2*(i*Width+j)+0]);
    printf("\n");
  }

  printf("Filtered v:\n");
  printf("\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.7f ",uFilt[2*(i*Width+j)+1]);
    printf("\n");
  }

  {
    double sum;
    sum=0.;
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1)
        sum += uFilt[2*(i*Width+j)+0];
    printf("sum=%.22f",sum);
  }

  free(uField);
  free(uBuff);
  free(uFilt);
  free(X);
  free(Y);
  free(Xbuff);
  free(Ybuff);
  
  return 0;
}