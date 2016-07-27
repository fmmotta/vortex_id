#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "floodFill.h"
#include "lambdaInit.h"
#include "stencilExtended.h"

#define filtered true

#define dbgPrint(num,num2) if(DEBUG_PRINT) printf("check point - %d-%d\n",(num),(num2))

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));\
                                  if(ptr==NULL){                         \
                                    printf("memory not allocked\n");     \
                                    return 1;                            \
                                  }                                      \
                                  else{                                  \
                                    for(i=0;i<(size);i+=1)               \
                                      ptr[i]=(type) 0;                   \
                                  }                                      \

int main(int argc,char **argv){
  int Width = 200, Height = 200;
  int i,j,err, padWidth=2,nVortex=3;
  int type = 1; // 0 = swst, 1 = vc
  int label[Width*Height],**eqClass;
  double *parVortex,x0[2],dx[2],xf[2],*sField=NULL;
  double *gField=NULL,*g2Field,*uField=NULL,X[Width],Y[Height];
  double *ux,*uy,*uxxx,*uyyy,*uxxy,*uxyy;
  double *uBuff=NULL,Xbuff[Width+4],Ybuff[Height+4];
  double x,y,v0y0 = 0.00;
  FILE *vFile;

  if(argc!=3){
    printf("Need number of vortices and file with vortices\n");
    return 1;
  }

  nVortex = atoi(argv[1]);
  vFile = fopen(argv[2],"r");
  if(vFile==NULL){
    printf("could not open vortex file\n");
    return 2;
  }

  fieldAlloc(sField,Height*Width,double);
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(ux,2*Height*Width,double);
  fieldAlloc(uy,2*Height*Width,double);
  fieldAlloc(uxxx,2*Height*Width,double);
  fieldAlloc(uxxy,2*Height*Width,double);
  fieldAlloc(uxyy,2*Height*Width,double);
  fieldAlloc(uyyy,2*Height*Width,double);
  fieldAlloc(uBuff,2*(Height+2*padWidth)*(Width+2*padWidth),double);
  fieldAlloc(gField,4*Height*Width,double);
  fieldAlloc(g2Field,4*Height*Width,double);
  fieldAlloc(parVortex,4*nVortex,double)

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }

  //x0[0]=-10.; xf[0]= 10.; 
  //x0[1]=-10.; xf[1]= 10.; 
  
  fscanf(vFile,"%d %d",&Height,&Width);
  fscanf(vFile,"%lf %lf",&(x0[0]),&(xf[0]));
  fscanf(vFile,"%lf %lf",&(x0[1]),&(xf[1]));
  
  dx[0] = (xf[0]-x0[0])/Width; dx[1] = (xf[1]-x0[1])/Height;

  for(i=0;i<nVortex;i+=1){
    fscanf(vFile,"%lf %lf %lf %lf", &(parVortex[4*i+0])
                                  , &(parVortex[4*i+1])
                                  , &(parVortex[4*i+2])
                                  , &(parVortex[4*i+3]));
  }
  fclose(vFile);
  /*
  parVortex[0+0]=1.; parVortex[0+1]=1.; parVortex[0+2]=-2.; parVortex[0+3]=0.;
  parVortex[4+0]=1.; parVortex[4+1]=1.;  parVortex[4+2]=2.; parVortex[4+3]=0.;
  parVortex[8+0]=1.; parVortex[8+1]=1.;  parVortex[8+2]=0.; parVortex[8+3]=4.;
  */
  for(j=0;j<Width;j+=1)
    X[j] = x0[0] + j*dx[0];

  for(i=0;i<Height;i+=1)
    Y[i] = x0[1] + i*dx[1];

  err = XtoXbuff(Width,X,Xbuff,2);
  err = XtoXbuff(Height,Y,Ybuff,2);
  
  err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
  if(err!=0)
    printf("Problems in addUSingleOseen\n");
  
  /* WARNING: come back here latter*/
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      uField[2*(i*Width+j)+0] += v0y0*Y[i];

  err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
  if(err!=0)
    printf("Problems in uFieldTouBuff\n");
  
  if(type == 0){
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
  }
  else{
    err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
    if(err!=0)
      printf("Problems in uFieldTouBuff\n");
  
    // \partial_x \vec u
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUx5point\n");

    // \partial_y \vec u
    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUy5point\n");

    // \partial_xxx \vec u
    err = UtoUxxx5point(Height,Width,uxxx,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUxxx5point\n");

    // \partial_yyy \vec u
    err = UtoUyyy5point(Height,Width,uyyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUyyy5point\n");
  
    // \partial_xxy \vec u
    err = uFieldTouBuff(Height,Width,uy,uBuff,padWidth);
    if(err!=0)
      printf("Problems in uFieldTouBuff\n");
  
    err = UtoUxx5point(Height,Width,uxxy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUxx5point\n");

    // \partial_yyx \vec u
    err = uFieldTouBuff(Height,Width,ux,uBuff,padWidth);
    if(err!=0)
      printf("Problems in uFieldTouBuff\n");

    err = UtoUyy5point(Height,Width,uxyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      printf("Problems in UtoUyy5point\n");

    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        gField[4*(i*Width+j)+0] = ux[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+1] = uy[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+2] = ux[2*(i*Width+j)+1];
        gField[4*(i*Width+j)+3] = uy[2*(i*Width+j)+1];
      }
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        g2Field[4*(i*Width+j)+0] = uxxy[2*(i*Width+j)+1]-uxyy[2*(i*Width+j)+0];
        g2Field[4*(i*Width+j)+1] = uxyy[2*(i*Width+j)+1]-uyyy[2*(i*Width+j)+0];
        g2Field[4*(i*Width+j)+2] = uxxy[2*(i*Width+j)+0]-uxxx[2*(i*Width+j)+1];
        g2Field[4*(i*Width+j)+3] = uxyy[2*(i*Width+j)+0]-uxxy[2*(i*Width+j)+1];
      }
    
    if(filtered)
      err=gradU2UtoLambda(Height,Width,gField,g2Field,&sField);
    else
      err=gradUtoLamb(Height,Width,g2Field,&sField);
    if(err!=0)
      printf("Problems in gradU2UtoLambda\n");
  }

  err = floodFill(sField,Width,Height,eqClass,label);
  if(err!=0)
    printf("Problems in floodFill\n");

  err = renameLabels(Width,Height,label);
  if(err>0)
    printf("%d connected component(s)\n",err);
  else
    printf("problems with renameLabels\n");

  {
    FILE *dadosout;
    dadosout=fopen("data/sField.dat","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        x = x0[0] + j*dx[0];
        y = x0[1] + i*dx[1];
        
        fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
      }
    fclose(dadosout);dadosout=NULL;

    dadosout=fopen("data/labelField.dat","w");
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

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  return 0;
}
