#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"
#include "lambdaInit.h"
#include "stencilManip.h"

int main(int argc,char **argv){
  const int Width = 100, Height = 100, Pop=10,nVortex=3;
  int i,j,err,ngbr,found,type=0;
  int nbList[8],label[Width*Height],eqList[Pop],**eqClass;
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float *gField=NULL,*uField=NULL,X[Width],Y[Height];
  float x,y,v0y0 = 0.0;

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  x0[0]=-5.; xf[0]= 5.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=-5.; xf[1]= 5.; dx[1] = (xf[1]-x0[1])/Width;

  parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=-2.; parVortex[3]=0.;
  parVortex[4+0]=1.; parVortex[4+1]=1.; parVortex[4+2]=2.; parVortex[4+3]=0.;
  parVortex[8+0]=1.; parVortex[8+1]=1.; parVortex[8+2]=0.; parVortex[8+3]=4.;
 
  gField = (float *)malloc(4*Height*Width*sizeof(float));
  if(gField==NULL){
    printf("memory not allocked\n");
    return 1;
  }
  for(i=0;i<4*Height*Width;i+=1)
    gField[i]=0.;

  uField = (float *)malloc(2*Height*Width*sizeof(float));
  if(uField==NULL){
    printf("memory not allocked\n");
    return 1;
  }
  for(i=0;i<2*Height*Width;i+=1)
    uField[i]=0.;

  for(j=0;j<Width;j+=1)
    X[j] = x0[0] + j*dx[0];

  for(i=0;i<Height;i+=1)
    Y[i] = x0[1] + i*dx[1];
  /*
  err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&g2Field);
  if(err!=0)
    printf("Problems in addSingleOseen\n");*/

  err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
  if(err!=0)
    printf("Problems in addUSingleOseen\n");
  
  err = UToGradUnUnifFullFrame(Height,Width,type,x0,dx,X,Y,uField,gField);
  if(err!=0)
    printf("Problems in uFieldToGradUopenFOAM\n");

  err = gradUtoLamb(Height,Width,gField,&sField);
  if(err!=0)
    printf("Problems in gradUtoLamb\n");
  
  err = floodFill(sField,Width,Height,eqClass,label);
  if(err!=0)
    printf("Problems in floodFill\n");

  err = renameLabels(Height,Width,label);
  if(err>0)
    printf("%d connected component(s)\n",err);
  else
    printf("problems with renameLabels - %d\n",err);

  {
    FILE *dadosout;
    dadosout=fopen("data/initULambOseen2D-4.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
      }

    fclose(dadosout);dadosout=NULL;

    dadosout=fopen("data/labelULambOseen2D-4.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %2d \n",x,y,label[i*Width+j]+1);
      }
      fprintf(dadosout,"\n");
    }

    fclose(dadosout);

    dadosout=fopen("data/uULambOseen2D-4.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f \n",x,y,
                                          uField[2*(i*Width+j)+0],
                                          uField[2*(i*Width+j)+1]);
      }
      fprintf(dadosout,"\n");
    }
    fclose(dadosout);

    dadosout=fopen("data/gUField-4.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f %f %f \n",x,y,
                                          gField[4*(i*Width+j)+0],
                                          gField[4*(i*Width+j)+1],
                                          gField[4*(i*Width+j)+2],
                                          gField[4*(i*Width+j)+3]);
      }
    }
    fclose(dadosout);
  }
  
  if(sField!=NULL)
    free(sField);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  return 0;
}
