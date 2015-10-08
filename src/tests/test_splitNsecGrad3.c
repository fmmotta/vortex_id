#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"
#include "lambdaInit.h"
#include "stencilExtended.h"
#include "vortexExtraction.h"

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc(size*sizeof(type)); \
                                  if(ptr==NULL){                        \
                                    printf("memory not allocked\n");    \
                                    return 1;                           \
                                  }                                     \
                                  else{                                 \
                                    for(i=0;i<size;i+=1)                \
                                      ptr[i]=(type) 0;                  \
                                  }                                     \

int oseenUxxy(int nVortex,double *parVortex, double *x0, double *dx, 
              int Height,int Width, double **sRefOut);

int main(int argc,char **argv){
  const int Width = 100, Height = 100, Pop=10,nVortex=3;
  int i,j,err,ngbr,found, padWidth=2;
  int nbList[8],label[Width*Height],eqList[Pop],**eqClass;
  double parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  double *gField=NULL,*g2Field=NULL,*uField=NULL,X[Width],Y[Height];
  double *uBuff=NULL,Xbuff[Width+4],Ybuff[Height+4],*g2Ref;
  double *ux,*uy,*uxxy,*uxyy,*uxxx,*uyyy,*w,*sRef1,*sRef2;
  double x,y,v0y0 = 0.0,theta;

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

  //parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=0.; parVortex[3]=0.;
  parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=-2.; parVortex[3]=0.;
  parVortex[4+0]=1.; parVortex[4+1]=1.; parVortex[4+2]=2.; parVortex[4+3]=0.;
  parVortex[8+0]=1.; parVortex[8+1]=1.; parVortex[8+2]=0.; parVortex[8+3]=4.;
 
  fieldAlloc(sField ,Height*Width,double);
  fieldAlloc(sRef1 ,Height*Width,double);
  fieldAlloc(sRef2 ,Height*Width,double);
  fieldAlloc(gField ,4*Height*Width,double);
  fieldAlloc(g2Field,4*Height*Width,double);
  fieldAlloc(g2Ref,4*Height*Width,double);
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(  ux  ,2*Height*Width,double);
  fieldAlloc(  uy  ,2*Height*Width,double);
  fieldAlloc( uxxy ,2*Height*Width,double);
  fieldAlloc( uxyy ,2*Height*Width,double);
  fieldAlloc( uxxx ,2*Height*Width,double);
  fieldAlloc( uyyy ,2*Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);

  for(j=0;j<Width;j+=1)
    X[j] = x0[0] + ((double)j)*dx[0];
  for(i=0;i<Height;i+=1)
    Y[i] = x0[1] + ((double)i)*dx[1];

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
  if(err!=0)
    printf("Problems in addUSingleOseen\n");
  
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
    printf("Problems in UtoUx5point\n");

  // \partial_yyy \vec u
  err = UtoUyyy5point(Height,Width,uyyy,uBuff,Xbuff,Ybuff);
  if(err!=0)
    printf("Problems in UtoUy5point\n");
  
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

  //err = gradUtoLamb(Height,Width,gField,&sField);
  err = gradUtoLamb(Height,Width,g2Field,&sField);
  //err=gradU2UtoLambda(Height,Width,gField,g2Field,&sField);
  if(err!=0)
    printf("Problems in gradU2UtoLambda\n");

  //theta = 0.; // Seems to be the minimal threshold to recover the correct
  //err=applySwirlingStrengthThreshold(Height,Width,sField,theta);// behaviour
  //if(err!=0)
  //  printf("Problems in applySwirlingStrengthThreshold\n");

  err = floodFill(sField,Width,Height,eqClass,label);
  if(err!=0)
    printf("Problems in floodFill\n");

  err = renameLabels(Height,Width,label);
  if(err>0)
    printf("%d connected component(s)\n",err);
  else
    printf("problems with renameLabels - %d\n",err);

  err = addOseen2ndGrad(nVortex,parVortex,x0,dx,Height,Width,&g2Ref);
  if(err!=0)
    printf("problems calculating g2Ref\n");

  err = gradUtoLamb(Height,Width,g2Ref,&sRef2);
  //err=s2ndGradUtoLamb(nVortex,parVortex,x0,dx,Height,Width,g2Ref,sRef2);
  if(err!=0)
    printf("Problems in gradU2UtoLambda\n");

  {
    FILE *dadosout;
    dadosout=fopen("data/initUSplit-3.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %.12f %.12f\n",x,y,sField[i*Width+j]
                                              ,sRef2[i*Width+j]);
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
    
    err = oseenUxxy(nVortex,parVortex,x0,dx,Height,Width,&sRef1);
    if(err!=0)
      printf("problems calculating ossen Uxxy\n");
    dadosout=fopen("data/uxxyUsplit-3.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f ",x,y,
                                          uxxy[2*(i*Width+j)+0],
                                          uxxy[2*(i*Width+j)+1]);
        fprintf(dadosout,"%f\n",sRef1[i*Width+j]);
      }
      fprintf(dadosout,"\n");
    }
    fclose(dadosout);

    dadosout=fopen("data/uxyyUsplit-3.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f \n",x,y,
                                          uxyy[2*(i*Width+j)+0],
                                          uxyy[2*(i*Width+j)+1]);
      }
      fprintf(dadosout,"\n");
    }
    fclose(dadosout);

    dadosout=fopen("data/uxxxUsplit-3.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f \n",x,y,
                                          uxxx[2*(i*Width+j)+0],
                                          uxxx[2*(i*Width+j)+1]);
      }
      fprintf(dadosout,"\n");
    }
    fclose(dadosout);

    dadosout=fopen("data/uyyyUsplit-3.txt","w");
    for(i=0;i<Height;i+=1){
      for(j=0;j<Width;j+=1){
        y = x0[0] + i*dx[0];
        x = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f %f \n",x,y,
                                          uyyy[2*(i*Width+j)+0],
                                          uyyy[2*(i*Width+j)+1]);
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
        
        fprintf(dadosout,"%f %f %.12f %.12f %.12f %.12f ",x,y,
                                          g2Field[4*(i*Width+j)+0],
                                          g2Field[4*(i*Width+j)+1],
                                          g2Field[4*(i*Width+j)+2],
                                          g2Field[4*(i*Width+j)+3]);
        fprintf(dadosout,"%.12f %.12f %.12f %.12f \n",g2Ref[4*(i*Width+j)+0],
                                          g2Ref[4*(i*Width+j)+1],
                                          g2Ref[4*(i*Width+j)+2],
                                          g2Ref[4*(i*Width+j)+3]);
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

int oseenUxxy(int nVortex,double *parVortex, double *x0, double *dx, 
              int Height,int Width, double **sRefOut){
  int i,j,k;
  double *sRef;
  double s,a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  double s0,xa2,yb2,R2;

  if(*sRefOut==NULL)
    return 1;
  sRef = *sRefOut;
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      s=0.;

      y = x0[0] + i*dx[0]; 
      x = x0[1] + j*dx[1]; 
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        xa2=(x-a)*(x-a);
        yb2=(y-b)*(y-b);
        R2 = R*R;
        
        // added a clause for small r/R
        // future : review gradU designation
        if(r<=0)
          s+=G/(2.*M_PI*R2*R2);
        else if((r>0)&&(r/R<cutoff))
          s+=(G/(2.*M_PI*R2*R2))*(1.-2.*(r2/R2));
        else{          
          fa = -(G*exp(-r2/R2))/(M_PI*R2*R2*R2*r2*r2);
          fb = 3.*(exp(r2/R2)-1.)*R2*R2*R2;
          
          s0  = 4.*xa2*yb2*r2*r2*r2;
          s0 += fb*(xa2*xa2 - 6.*xa2*yb2 + yb2*yb2);
          s0 += 2.*R2*r2*r2*(xa2*xa2 - 4.*xa2*yb2 + yb2*yb2);
          s0 += -3.*R2*R2*(xa2*xa2*xa2-5.*xa2*xa2*yb2-5.*xa2*yb2*yb2+yb2*yb2*yb2);
          
          s += s0*fa;
        }
      }

      sRef[i*Width+j]+=s;
    }

  return 0;
}