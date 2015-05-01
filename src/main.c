#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"
#include "lambdaInit.h"

int main(int argc,char **argv){
  const int Width = 100, Height = 100, Pop=10,nVortex=3;
  int i,j,err,ngbr,found;
  int nbList[8],label[Width*Height],eqList[Pop];
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float x,y,v0y0 = 0.05;
  
  /*               
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
  }*/
  

  //x0[0]=-4.75; xf[0]= 5.25; dx[0] = (xf[0]-x0[0])/Height;
  //x0[1]=-4.75; xf[1]= 5.25; dx[1] = (xf[1]-x0[1])/Width;
  //x0[0]=-5.; xf[0]=5.; dx[0] = (xf[0]-x0[0])/Height;
  //x0[1]= 0.; xf[1]=0.; dx[1] = 0.;
  
  x0[0]=-5.; xf[0]= 5.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=-5.; xf[1]= 5.; dx[1] = (xf[1]-x0[1])/Width;

  parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=-2.; parVortex[3]=0.;
  parVortex[4+0]=1.; parVortex[4+1]=1.; parVortex[4+2]=2.; parVortex[4+3]=0.;
  parVortex[8+0]=1.; parVortex[8+1]=1.; parVortex[8+2]=0.; parVortex[8+3]=4.;
 
  //parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=0.; parVortex[3]=0.;
  
  printf("initOseenShear2D\n");
  err = initOseenShear2D(nVortex,parVortex,x0,dx,Height,Width,v0y0,&sField);
  
  printf("floodFill\n");
  err = floodFill(sField,Width,Height,label);

  {
    FILE *dadosout;
    dadosout=fopen("data/initLambOseen2D.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        x = x0[0] + j*dx[0];
        y = x0[1] + i*dx[1];
        
        fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
      }

    fclose(dadosout);dadosout=NULL;

    dadosout=fopen("data/labelLambOseen2D.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        x = x0[0] + j*dx[0];
        y = x0[1] + i*dx[1];
        
        fprintf(dadosout,"%f %f %2d \n",x,y,label[i*Width+j]+1);
      }
      fprintf(dadosout,"\n");
    }

    fclose(dadosout);
  }
  
  if(sField!=NULL)
    free(sField);
  return 0;
}
