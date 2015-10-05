#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"
#include "lambdaInit.h"
#include "stencilExtended.h"

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc(size*sizeof(type));if(ptr==NULL){printf("memory not allocked\n"); return 1;}

int main(int argc,char **argv){
  const int Width = 100, Height = 100, Pop=10,nVortex=3;
  int i,j,err,ngbr,found, padWidth=2;
  int nbList[8],label[Width*Height],eqList[Pop],**eqClass;
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float *gField=NULL,*uField=NULL,X[Width],Y[Height];
  float *uBuff=NULL,Xbuff[Width+4],Ybuff[Height+4];
  float *ux,*uy,*uxxy,*uxyy,*uxxx,*uyyy,*w;
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
 
  fieldAlloc(gField,4*Height*Width,float);
  for(i=0;i<4*Height*Width;i+=1)
    gField[i]=0.;

  fieldAlloc(uField,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uField[i]=0.;

  fieldAlloc(uBuff,2*(Height+2*padWidth)*(Width+2*padWidth),float);
  for(i=0;i<2*(Height+2*padWidth)*(Width+2*padWidth);i+=1)
    uBuff[i]=0;

  fieldAlloc(ux,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    ux[i]=0.;

  fieldAlloc(uy,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uy[i]=0.;

  fieldAlloc(uxxy,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uxxy[i]=0.;

  fieldAlloc(uxyy,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uxyy[i]=0.;

  fieldAlloc(uxxx,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uxxx[i]=0.;

  fieldAlloc(uyyy,2*Height*Width,float);
  for(i=0;i<2*Height*Width;i+=1)
    uyyy[i]=0.;

  for(j=0;j<Width;j+=1)
    X[j] = x0[0] + j*dx[0];
  for(i=0;i<Height;i+=1)
    Y[i] = x0[1] + i*dx[1];

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  /*
  err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&g2Field);
  if(err!=0)
    printf("Problems in addSingleOseen\n");*/

  err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
  if(err!=0)
    printf("Problems in addUSingleOseen\n");
  
  err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
    printf("Problems in uFieldTouBuff\n");
  
  /*
  err = UToGradUnUnifFullFrame(Height,Width,type,x0,dx,X,Y,uField,gField);
  if(err!=0)
    printf("Problems in uFieldToGradUopenFOAM\n");*/

  err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
  if(err!=0)
    printf("Problems in UtoUx5point\n");

  err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
  if(err!=0)
    printf("Problems in UtoUy5point\n");

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gField[4*(i*Width+j)+0] = ux[2*(i*Width+j)+0];
      gField[4*(i*Width+j)+1] = uy[2*(i*Width+j)+0];
      gField[4*(i*Width+j)+2] = ux[2*(i*Width+j)+1];
      gField[4*(i*Width+j)+3] = uy[2*(i*Width+j)+1];
    }

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
    dadosout=fopen("data/initUSplit-3.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
      }

    fclose(dadosout);dadosout=NULL;

    dadosout=fopen("data/labelUSplit-3.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %2d \n",x,y,label[i*Width+j]+1);
      }
      fprintf(dadosout,"\n");
    }

    fclose(dadosout);

    dadosout=fopen("data/uUsplit-3.txt","w");
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

    dadosout=fopen("data/uUbuffsplit-3.txt","w");
    for(i=0;i<(Height+2*padWidth);i+=1){
      for(j=0;j<(Width+2*padWidth);j+=1){
        y = Ybuff[i];
        x = Xbuff[j];
        
        fprintf(dadosout,"%f %f %f %f \n",x,y,
                                          uBuff[2*(i*(Width+2*padWidth)+j)+0],
                                          uBuff[2*(i*(Width+2*padWidth)+j)+1]);
      }
      fprintf(dadosout,"\n");
    }
    fclose(dadosout);

    dadosout=fopen("data/gUSplit-3.txt","w");
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
