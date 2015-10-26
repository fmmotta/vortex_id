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
#include "inputManager.h"
#include "essayHandler.h"

int fprintfRunParamSigned(FILE *dadosgen,long long int seed,double x0[],
                         double xf[],double dx[],double Gmin, double Gmax,
                         double rmin, double rmax, double xmin[], double xmax[], 
                         double v0y0)
{
  int i;

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
                   gsl_histogram *ia, gsl_histogram *ib){
  int i;

  for(i=0;i<nVortex;i+=1){
    gsl_histogram_increment(iG,parVortex[4*i+0]);
    gsl_histogram_increment(iRc,parVortex[4*i+1]);
    gsl_histogram_increment(ia,parVortex[4*i+2]);
    gsl_histogram_increment(ib,parVortex[4*i+3]);
  }

  return 0;
}

int fprintVortex(FILE *dadosout, int run,int nVortex, double *vCatalog){
  int i;

  if(dadosout==NULL || run<0 || nVortex<=0 || vCatalog==NULL)
    return 1;

  for(i=0;i<nVortex;i+=1)
    fprintf(dadosout,"%d %d : %f %f %f %f\n",run,i,vCatalog[4*i+0]
                                                  ,vCatalog[4*i+1]
                                                  ,vCatalog[4*i+2]
                                                  ,vCatalog[4*i+3]);

  fprintf(dadosout,"\n");

  return 0;
}

int fprintsField(FILE *dadosout,double *x0,double *dx,
                 int Width, int Height, double *sField){
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

  if(runType==0){
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
  else if(runType==1){
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
  else if(runType==2){
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
  else{
    printf("Non-Identified run-type - %d\n",runType);
    return -2;
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