#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_histogram.h>
#include "ini.h"
#include "mt64.h"
#include "vortexGen.h"
#include "floodFill.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "stencilExtended.h"
#include "vortexExtractionExtend.h"
#include "inputManager.h"
#include "essayHandler.h"

#define filtered true

int fprintfRunParamSigned(FILE *dadosgen,long long int seed,double x0[],
                         double xf[],double dx[],double Gmin, double Gmax,
                         double rmin, double rmax, double xmin[], double xmax[], 
                         double v0y0)
{
  if(dadosgen==NULL)
    return 1;

  fprintf(dadosgen,"seed: %lld\n",seed);
  fprintf(dadosgen,"\ndomain: xi xf dx\n");
  fprintf(dadosgen,"%f %f %f\n",x0[0],xf[0],dx[0]);
  fprintf(dadosgen,"%f %f %f\n",x0[1],xf[1],dx[1]);
  fprintf(dadosgen,"\nvortex params: Uniform\n");
  fprintf(dadosgen,"G :+- %f %f\n",Gmin,Gmax);
  fprintf(dadosgen,"Rc: %f %f\n",rmin,rmax);
  fprintf(dadosgen,"a : %f %f\n",xmin[0],xmax[0]);
  fprintf(dadosgen,"b : %f %f\n",xmin[1],xmax[1]);
  fprintf(dadosgen,"\nshear v0/y0=%f\n",v0y0);
  
  return 0;
}

int histoIncVortex(int nVortex, double *parVortex,
                   gsl_histogram *iG, gsl_histogram *iRc,
                   gsl_histogram *ia, gsl_histogram *ib)
{
  int i;

  for(i=0;i<nVortex;i+=1){
    gsl_histogram_increment(iG,parVortex[4*i+0]);
    gsl_histogram_increment(iRc,parVortex[4*i+1]);
    gsl_histogram_increment(ia,parVortex[4*i+2]);
    gsl_histogram_increment(ib,parVortex[4*i+3]);
  }

  return 0;
}

int fprintVortex(FILE *dadosout, int run,int nVortex, double *vCatalog)
{
  int i;

  if(dadosout==NULL || run<0 || nVortex<=0 || vCatalog==NULL)
    return 1;

  for(i=0;i<nVortex;i+=1)
    fprintf(dadosout,"%f %f %f %f\n",vCatalog[4*i+0]
                                    ,vCatalog[4*i+1]
                                    ,vCatalog[4*i+2]
                                    ,vCatalog[4*i+3]);

  fprintf(dadosout,"\n");

  return 0;
}

int fprintSafeVortex(FILE *dadosout, int run,int nVortex, 
                     double *vCatalog,int Height,int Width,
                     double *X,double *Y)
{
  int i;

  if(dadosout==NULL || run<0 || nVortex<=0 || vCatalog==NULL){
    printf("input: %p %d %d %p\n",dadosout,run,nVortex,vCatalog);
    return 1;
  }

  for(i=0;i<nVortex;i+=1){
    if(vCatalog[12*i+2] < X[0] || vCatalog[12*i+2] > X[Width-1] ||
       vCatalog[12*i+3] < Y[0] || vCatalog[12*i+3] > Y[Height-1])
      continue;
    fprintf(dadosout,"%f %f %f %f\n",vCatalog[12*i+0]
                                    ,vCatalog[12*i+1]
                                    ,vCatalog[12*i+2]
                                    ,vCatalog[12*i+3]);
  }
  fprintf(dadosout,"\n");

  return 0;
}

int fprintSafeVortexMoments(FILE *dadosout, int run,int size,int nVortex, 
                            double *rCatalog,int Height,int Width,
                            double *X,double *Y)
{
  int i;
  double Xcm,Ycm,rg2,a,b,c,d,norm,ac,bc;
  double L1,L2,trM2,detM,gradU[2][2];

  if(dadosout==NULL || run<0 || nVortex<=0 || rCatalog==NULL){
    printf("input: %p %d %d %p\n",dadosout,run,nVortex,rCatalog);
    return 1;
  }

  for(i=0;i<nVortex;i+=1){
    if(rCatalog[size*i+2] < X[0] || rCatalog[size*i+2] > X[Width-1] ||
       rCatalog[size*i+3] < Y[0] || rCatalog[size*i+3] > Y[Height-1])
      continue;
    Xcm = rCatalog[size*i+2];
    Ycm = rCatalog[size*i+3];

    // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html

    a = rCatalog[size*i+4] - Xcm*Xcm;
    b = rCatalog[size*i+5] - Xcm*Ycm;
    c = rCatalog[size*i+6] - Ycm*Xcm;
    d = rCatalog[size*i+7] - Ycm*Ycm;

    fprintf(dadosout,"%f %f %f %f ",rCatalog[size*i+0]
                                   ,rCatalog[size*i+1]
                                   ,rCatalog[size*i+2]
                                   ,rCatalog[size*i+3]);

    fprintf(dadosout,"%.8f %.8f %.8f %.8f ",a,b,c,d);

    rg2 = a+d;
    fprintf(dadosout,"%f ", sqrt(rg2));
    
    trM2  = (a+d)/2.;
    detM  = a*d-b*c;

    L1 = trM2+sqrt(trM2*trM2-detM);
    L2 = trM2-sqrt(trM2*trM2-detM);
    
    fprintf(dadosout,"%f %f ",L1,L2);

    if(a>0.){
      
      if(b>0){
        norm = sqrt(b*b+(L1-a)*(L1-a));
        fprintf(dadosout,"%f %f ",b/norm,(L1-a)/norm);
       
        norm = sqrt(b*b+(L2-a)*(L2-a));
        fprintf(dadosout,"%f %f ",b/norm,(L2-a)/norm);
      }
      else{
        norm = sqrt(b*b+(L1-a)*(L1-a));
        fprintf(dadosout,"%f %f ",-b/norm,-(L1-a)/norm);
       
        norm = sqrt(b*b+(L2-a)*(L2-a));
        fprintf(dadosout,"%f %f ",-b/norm,-(L2-a)/norm); 
      }
    }
    else{
      fprintf(dadosout,"%f %f ",1.,0.);
      fprintf(dadosout,"%f %f ",0.,1.);
    }
    
    norm = sqrt(b*b+(L1-a)*(L1-a));
    if(b>0)
      fprintf(dadosout,"%f ",(L1-a)/norm);
    else
      fprintf(dadosout,"%f ",-(L1-a)/norm);
    
    gradU[0][0] = rCatalog[size*i+8];
    gradU[0][1] = rCatalog[size*i+9];
    gradU[1][0] = rCatalog[size*i+10];
    gradU[1][1] = rCatalog[size*i+11];

    fprintf(dadosout,"%f %f %f %f ",gradU[0][0],
                                    (gradU[0][1]+gradU[1][0])/2.,
                                    (gradU[0][1]+gradU[1][0])/2.,
                                    gradU[1][1]);
    
    ac = (gradU[0][1]+gradU[1][0])/2.;
    bc = gradU[1][1];
    fprintf(dadosout,"%f %f",gradU[1][0]-gradU[0][1],sqrt(ac*ac+bc*bc));
    
    fprintf(dadosout,"\n");
  }


  return 0;
}

int fprintsField(FILE *dadosout,double *x0,double *dx,
                 int Width, int Height, double *sField)
{
  int i,j;
  double x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || sField==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      y = x0[0] + i*dx[0];
      x = x0[1] + j*dx[1];

      fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
    }

  return 0;
}


int fprintUsfield(FILE *dadosout,double *X,double *Y, 
                  int Height,int Width, double *sField)
{
  int i,j;
  double x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || sField==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      y = Y[i];
      x = X[j];

      //fprintf(dadosout,"%f %f %.12f\n",x,y,sField[i*Width+j]);
      fprintf(dadosout,"%f %f %.12f\n",x,y,log(1.+sField[i*Width+j]));
    }

  return 0;
}

int fprintLabels(FILE *dadosout,double *x0,double *dx,
                 int Width, int Height, int *label){
  int i,j;
  double x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || label==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      y = x0[0] + i*dx[0];
      x = x0[1] + j*dx[1];

      fprintf(dadosout,"%f %f %d\n",x,y,label[i*Width+j]);
    }  

  return 0;
}

int fprintUlabels(FILE *dadosout,double *X,double *Y,
                  int Height, int Width,int *label)
{
  int i,j;
  double x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || label==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      y = Y[i];
      x = X[j];

      fprintf(dadosout,"%f %f %d\n",x,y,label[i*Width+j]);
    }  

  return 0;
}

int fprintUpresence(FILE *dadosout,double *X,double *Y,
                  int Height, int Width,int *label)
{
  int i,j,p;
  double x,y;

  if(dadosout==NULL || Width<0 || Height<=0 || label==NULL)
    return 1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      y = Y[i];
      x = X[j];
      
      if(label[i*Width+j]>=0)
        p=1;
      else
        p=0;
      fprintf(dadosout,"%f %f %d\n",x,y,p);
    }  

  return 0;
}

//---------------------------------------------------------------------

void printScalarFields(int Height,int Width, double *X,double *Y,int n,
                     char *folder,int index,double *sField,int *label)
{
  char filename[512+1];
  FILE *dadosout;

  sprintf(filename,"%s/sField-p%d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUsfield(dadosout,X,Y,Height,Width,sField);
  fclose(dadosout);

  sprintf(filename,"%s/labels-p%d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUlabels(dadosout,X,Y,Height,Width,label);
  fclose(dadosout);

  sprintf(filename,"%s/presence-p%d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUpresence(dadosout,X,Y,Height,Width,label);
  fclose(dadosout);
}

//----------------------------------------------------------------------

void printVorticesAndMoments(int Height,int Width, double *X,double *Y,double t,int n,
                             char *folder,int index,int nCnect,double *vCatalog,double *rCatalog)
{
  int err,dataSize=12;
  char filename[512+1];
  FILE *dadosVout;

  sprintf(filename,"%s/vortices-p%d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");   
  err=fprintVortex(dadosVout,n,nCnect,vCatalog);
  if(err!=0){printf("problems vortices\n"); exit(-6);}
  fclose(dadosVout);
      
  sprintf(filename,"%s/vorticesSafe-p%d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");
  err=fprintSafeVortex(dadosVout,n,nCnect,vCatalog,Height,Width,X,Y);
  if(err!=0){printf("problems vorticesSafe\n"); exit(-6);}
  fclose(dadosVout);
      
  sprintf(filename,"%s/vortexSafeMoments-p%d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");
  err=fprintSafeVortexMoments(dadosVout,n,dataSize,nCnect,rCatalog,Height,Width,X,Y);
  if(err!=0){printf("problems vortexSafeMoments\n"); exit(-6);}
  fclose(dadosVout);
}

//-----------------------------------------------------------------------

int genVortices(int genType,long long int seed, double xmin[],double xmax[], 
                int nFixVortex, double **parVortex,
                double Gmin,double Gmax,double rmin,double rmax,
                double numG,double numRc, double *Glist,double *Rclist)
{
  int err,nVortex;

  if(genType==0){
    nVortex = nFixVortex;
    err=genLOseenUniformList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                             nVortex,parVortex);
    if(err<0)
      return err;
    else if((err>0) && (err<nVortex))
      nVortex = err;
  }
  else if(genType==1){
    nVortex = nFixVortex;
    err=genLOseenNaryList(numG,Glist,numRc,Rclist,xmin,xmax,
                          seed,nVortex,parVortex);
    if(err<0)
      return err;
    else if((err>0) && (err<nVortex))
      nVortex = err;
  }
  else if(genType==2){
    nVortex = nFixVortex;
    err=genLOseenSignUniformList(Gmin,Gmax,rmin,rmax,xmin,xmax,seed,
                                 nVortex,parVortex);
    if(err<0)
      return err;
    else if((err>0) && (err<nVortex))
      nVortex = err;

    return nVortex;
  }
  else{
    printf("Non-Identified vortex generation type\n");
    return -1;
  }

  return 0;
}

int calcScalarField(int runType,int Height,int Width,double x0[],double dx[],
                    int nVortex,double *parVortex,double *gField,double v0y0,
                    double *sField)
{
  int err;

  if(runType==2){
    err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0){
      printf("Problems in addSingleOseen\n");
      return err;
    }

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0){
      printf("Problems in gradUtoLamb\n");
      return err;
    }
  }
  else if(runType==3){
    err = addOseen2ndGrad(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0){
      printf("Problems in addSingleOseen\n");
      return err;
    }

    err = s2ndGradUtoLamb(nVortex,parVortex,x0,dx,Height,Width,gField,sField);
    if(err!=0){
      printf("Problems in gradUtoLamb\n");
      return err;
    }
  }
  else if(runType==0){
    err = addSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0){
      printf("Problems in addSingleOseen\n");
      return err;
    }

    err=addConstXYShear(x0,dx,Height,Width,v0y0,&gField);
    if(err!=0)
      printf("Problems in addConstXYShear\n");

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0){
      printf("Problems in gradUtoLamb\n");
      return err;
    }
  }
  else if(runType==1){
    err = addOseen2ndGrad(nVortex,parVortex,x0,dx,Height,Width,&gField);
    if(err!=0){
      printf("Problems in addSingleOseen\n");
      return err;
    }

    err = s2ndGradUtoLambShear(nVortex,parVortex,x0,dx,Height,Width,
                               gField,sField,v0y0);
    if(err!=0){
      printf("Problems in gradUtoLamb\n");
      return err;
    }
  }
  else{
    printf("Non-Identified run-type - %d\n",runType);
    return -2;
  }
  return 0;
}

int calcUScalarField(int runType,int Height,int Width,int padWidth, 
                     double x0[],double dx[],double *X,double *Y, double *Xbuff,
                     double *Ybuff,int nVortex,double *parVortex, 
                     double *uField, double *uBuff,double *ux,double *uy,
                     double *uxxx, double *uyyy,double *uxxy, double *uxyy,
                     double *gField,double *g2Field,double v0y0,double *sField)
{
  int i,j,err;

  if(runType==0){
    err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
    if(err!=0)
      return -1;
    
    /* WARNING: come back here latter*/
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1)
        uField[2*(i*Width+j)+0] += v0y0*Y[i];
  
    err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
    if(err!=0)
      return -2;
  
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -3;

    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -4;

    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        gField[4*(i*Width+j)+0] = ux[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+1] = uy[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+2] = ux[2*(i*Width+j)+1];
        gField[4*(i*Width+j)+3] = uy[2*(i*Width+j)+1];
      }

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0)
      return -5;
  }
  else if(runType==1){
    
    err = addUSingleOseen(nVortex,parVortex,x0,dx,Height,Width,&uField);
    if(err!=0)
      return -1;
    
    /* WARNING: come back here latter*/
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1)
        uField[2*(i*Width+j)+0] += v0y0*Y[i];
      
    err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
    if(err!=0)
      return -2;
  
    // \partial_x \vec u
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -3;

    // \partial_y \vec u
    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -4;

    // \partial_xxx \vec u
    err = UtoUxxx5point(Height,Width,uxxx,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -5;

    // \partial_yyy \vec u
    err = UtoUyyy5point(Height,Width,uyyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -6;
  
    // \partial_xxy \vec u
    err = uFieldTouBuff(Height,Width,uy,uBuff,padWidth);
    if(err!=0)
      return -7;
  
    err = UtoUxx5point(Height,Width,uxxy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -8;

    // \partial_yyx \vec u
    err = uFieldTouBuff(Height,Width,ux,uBuff,padWidth);
    if(err!=0)
      return -9;

    err = UtoUyy5point(Height,Width,uxyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -10;

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
      return -11;
  }
  else{
    printf("Non-Identified run-type - %d\n",runType);
    return -20;
  }

  return 0;
}

int vortexReconstruction(int runType,int Height, int Width, int nCnect, 
                          double x0[],double dx[],double *sField, 
                          double *gField,int *label,double **vCatalog)
{
  int err;

  if(runType==0 || runType==2){
    err=vortexExtraction(Height,Width,nCnect,x0,dx,sField,
                         gField,label,vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }
  }
  else if(runType==1){
    err=vortexExt2ndSwirl(Height,Width,nCnect,x0,dx,sField,
                         gField,label,vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }
  } 
  else{
    printf("Non-Identified run-type\n");
    return -2;
  }

  return 0;
}

int vortexUReconstruction(int runType,int Height, int Width, int nCnect, 
                         double *X,double *Y,double *sField, 
                         double *gField,int *label,double **vCatalog)
{
  int err;

  if(runType==0 || runType==2){
    err=vortexExtFromSwirlStr(Height,Width,nCnect,X,Y,sField,
                              gField,label,vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }
  }
  else if(runType==1){
    err=vortexExtFromVortCurv(Height,Width,nCnect,X,Y,sField,
                              gField,label,vCatalog);
    if(err!=0){
      printf("error on vortexExtraction - %d\n",err);
      return err; 
    }
  } 
  else{
    printf("Non-Identified run-type\n");
    return -2;
  }

  return 0;
}

int writeGnuplotScript(char *filename,char *folder,char *tag,
                       int nRuns,int nVortex){
  char gfile[200+1];
  sprintf(gfile,"%s/%s",folder,filename);
  FILE *dadosout = fopen(gfile,"w");
  if(dadosout==NULL)
    return -1;
  
  fprintf(dadosout,"set yr [0:]\n");
  fprintf(dadosout,"set yl 'Counting' \n");
  fprintf(dadosout,"set key top center\n");
  fprintf(dadosout,"set style fill solid border -1\n");

  fprintf(dadosout,"set xl '$\\Gamma$'\n");
  fprintf(dadosout,"set yr [0:]\n");
  fprintf(dadosout,"plot 'histoOuG-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoOuG-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$r_c$'\n");
  fprintf(dadosout,"plot 'histoOuRc-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoOuRc-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$a$'\n");
  fprintf(dadosout,"plot 'histoOua-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoOua-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$b$'\n");
  fprintf(dadosout,"plot 'histoOub-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoOub-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl 'N'\n");
  fprintf(dadosout,"plot 'histoOuN-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoOuN-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$\\Gamma$'\n");
  fprintf(dadosout,"plot 'histoInG-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoInG-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");
  
  fprintf(dadosout,"set xl '$r_c$'\n");
  fprintf(dadosout,"plot 'histoInRc-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoInRc-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$a$'\n");
  fprintf(dadosout,"plot 'histoIna-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoIna-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fprintf(dadosout,"set xl '$b$'\n");
  fprintf(dadosout,"plot 'histoInb-%s.txt' using "
                   "(($1+$2)/2):3 with boxes title "
                   "'%d events of %d vortices'\n",tag,nRuns,nVortex);
  fprintf(dadosout,"set term epslatex standalone color colortext 12\n");
  fprintf(dadosout,"set out 'histoInb-%s.tex'\n",tag);
  fprintf(dadosout,"replot\n");
  fprintf(dadosout,"set term qt\n");
  fprintf(dadosout,"replot\n");

  fclose(dadosout);

  sprintf(gfile,"%s/graph_script.sh",folder);
  dadosout = fopen(gfile,"w");
  fprintf(dadosout,"gnuplot %s\n",filename);
  fprintf(dadosout,"pdflatex histoOuG-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoOuRc-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoOua-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoOub-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoOuN-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoInG-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoInRc-%s.tex\n",tag);
  fprintf(dadosout,"pdflatex histoIna-%s.tex \n",tag);
  fprintf(dadosout,"pdflatex histoInb-%s.tex \n",tag);
  fclose(dadosout);

  return 0;
}


int foamScalarField(int runType,int Height,int Width,int padWidth, 
                    double *X,double *Y, double *Xbuff,double *Ybuff, 
                    double *uField, double *uBuff,double *ux,double *uy,
                    double *uxxx, double *uyyy,double *uxxy, double *uxyy,
                    double *gField,double *g2Field,double v0y0,double *sField)
{
  int i,j,err;

  if(runType==0){
  
    err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
    if(err!=0)
      return -2;
  
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -3;

    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -4;

    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        gField[4*(i*Width+j)+0] = ux[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+1] = uy[2*(i*Width+j)+0];
        gField[4*(i*Width+j)+2] = ux[2*(i*Width+j)+1];
        gField[4*(i*Width+j)+3] = uy[2*(i*Width+j)+1];
      }

    err = gradUtoLamb(Height,Width,gField,&sField);
    if(err!=0)
      return -5;
  }
  else if(runType==1){    
    err = uFieldTouBuff(Height,Width,uField,uBuff,padWidth);
    if(err!=0)
      return -2;
  
    // \partial_x \vec u
    err = UtoUx5point(Height,Width,ux,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -3;

    // \partial_y \vec u
    err = UtoUy5point(Height,Width,uy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -4;

    // \partial_xxx \vec u
    err = UtoUxxx5point(Height,Width,uxxx,uBuff,Xbuff,Ybuff);
    if(err!=0)    
      return -5;

    // \partial_yyy \vec u
    err = UtoUyyy5point(Height,Width,uyyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -6;
  
    // \partial_xxy \vec u
    err = uFieldTouBuff(Height,Width,uy,uBuff,padWidth);
    if(err!=0)
      return -7;
    err = UtoUxx5point(Height,Width,uxxy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -8;

    // \partial_yyx \vec u
    err = uFieldTouBuff(Height,Width,ux,uBuff,padWidth);
    if(err!=0)
      return -9;
    err = UtoUyy5point(Height,Width,uxyy,uBuff,Xbuff,Ybuff);
    if(err!=0)
      return -10;

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
      return -11;
  }
  else{
    printf("Non-Identified run-type - %d\n",runType);
    return -20;
  }

  return 0;
}