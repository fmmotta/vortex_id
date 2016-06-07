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
#include "stencilExtended.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "vortexExtractionExtend.h"
#include "inputManager.h"
#include "essayHandler.h"

#define DEBUG_MODE false
#define DEBUG_PRINT false
#define PRINT_MODE false

#define dbgPrint(num,num2) if(DEBUG_PRINT) printf("check point - %d-%d\n",(num),(num2))

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));\
                                  if(ptr==NULL){                         \
                                    printf("memory not allocked\n");     \
                                    return 1;                            \
                                  }                                      \
                                  else{                                  \
                                    for(i=0;i<(size);i+=1)               \
                                      ptr[i]=(type) 0;                   \
                                  }      

int main(int argc,char **argv){
  int err,n;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10;
  int Width = 100, Height = 100,nVortex=5,nFixVortex=5,nRuns=1000,count=0;
  double vG,vR,va,vb;
  double hGmin=0.,hGmax=0.,hRcMin=0.,hRcMax=0.;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2];
  char folder[100+1],tag[100+1],filename[400+1],bkgFile[400+1];
  gsl_histogram *hG,*hRc,*ha,*hb,*hN;
  gsl_histogram *iG,*iRc,*ia,*ib;
  FILE *totalVin,*totalVout,*dadosout;
  configVar cfg;
  
  /*********************************/

  dbgPrint(0,0);

  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  dbgPrint(1,0);

  /* Loading Configuration -- I need something more concise */

  err=initConfig(&cfg);

  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
  
  dbgPrint(2,0);

  if(DEBUG_MODE==true){
    err=printConfig(&cfg);
    if(err!=0)
      return err;
  }
  
  if(cfg.dim!=2){
    printf("Dimension is not 2 - Can't Follow\n");
    return 2;
  }
  
  dbgPrint(2,1);

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  if(cfg.bkgFile!=NULL)
    strcpy(bkgFile,cfg.bkgFile);
  else
    bkgFile[0]='\0';

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }
  
  dbgPrint(3,0);

  Width   = cfg.Width;
  Height  = cfg.Height;
  nRuns   = cfg.nRuns;
  nFixVortex = cfg.nVortex;
  
  dbgPrint(4,0);

  /************************************/

  x0[0]   = cfg.x0[0]; x0[1] = cfg.x0[1]; 
  xmin[0] = x0[0]+1; xmin[1] = x0[1]+1;
  
  xf[0]   = cfg.xf[0]; xf[1] = cfg.xf[1]; 
  xmax[0] = xf[0]-1; xmax[1] = xf[1]-1;
  
  dx[0] = (xf[0]-x0[0])/Height;
  dx[1] = (xf[1]-x0[1])/Width;
  
  /**********************************/

  dbgPrint(5,0);

  hNG     = cfg.hNG;
  hNRc    = cfg.hNRc;
  hNa     = cfg.hNa;
  hNb     = cfg.hNb;
  hNN     = cfg.hNN;
  hGmin   = cfg.hGmin;
  hGmax   = cfg.hGmax;
  hRcMin  = cfg.hRcMin;
  hRcMax  = cfg.hRcMax;
  
  dbgPrint(6,0);

  err=freeConfig(&cfg);if(err!=0) return err;

  dbgPrint(7,0);

  /* histogram preparation - begin */
  hG  = gsl_histogram_alloc(hNG);  gsl_histogram_set_ranges_uniform(hG,hGmin,hGmax);
  hRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(hRc,hRcMin,hRcMax);
  ha  = gsl_histogram_alloc(hNa);  gsl_histogram_set_ranges_uniform(ha,xmin[0],xmax[0]);
  hb  = gsl_histogram_alloc(hNb);  gsl_histogram_set_ranges_uniform(hb,xmin[1],xmax[1]);
  hN  = gsl_histogram_alloc(hNN);  gsl_histogram_set_ranges_uniform(hN,0,2*nVortex);

  iG  = gsl_histogram_alloc(hNG);  gsl_histogram_set_ranges_uniform(iG,hGmin,hGmax);
  iRc = gsl_histogram_alloc(hNRc); gsl_histogram_set_ranges_uniform(iRc,hRcMin,hRcMax);
  ia  = gsl_histogram_alloc(hNa);  gsl_histogram_set_ranges_uniform(ia,xmin[0],xmax[0]);
  ib  = gsl_histogram_alloc(hNb);  gsl_histogram_set_ranges_uniform(ib,xmin[1],xmax[1]);
  /* histogram preparation - end*/

  dbgPrint(8,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
  }

  dbgPrint(9,0);

  sprintf(filename,"%s/inputVortexes.txt",folder);
  totalVin = fopen(filename,"r");
  if(totalVin==NULL){printf("Problems to open input vortex list\n");}

  sprintf(filename,"%s/outputVortexes.txt",folder);
  totalVout = fopen(filename,"r");
  if(totalVout==NULL){printf("Problems to open output vortex list\n");}

  dbgPrint(10,0);
  
  /**********************************************************/
  
  n=0;
  do{
  	err=fscanf(totalVin,"%lf%lf%lf%lf",&vG,&vR,&va,&vb);
  	if(DEBUG_PRINT && n%10000==0){
      printf("%d vortices processed\n",n);
      printf("v: %lf %lf %lf %lf\n",vG,vR,va,vb);
    }

    n+=1;
    if(err!=4)
      continue;

    /* printing to histogram */

    gsl_histogram_increment(iG ,vG);
    gsl_histogram_increment(iRc,vR);
    gsl_histogram_increment(ia ,va);
    gsl_histogram_increment(ib ,vb);

    /* Preparing for printing */
  }while(!feof(totalVin));

  dbgPrint(10,1);
  
  count = 0;
  n=0;
  do{
    err=fscanf(totalVout,"%lf%lf%lf%lf",&vG,&vR,&va,&vb);
  	if(DEBUG_PRINT && n%10000==0){
      printf("%d vortices processed\n",n);
      printf("v: %lf %lf %lf %lf\n",vG,vR,va,vb);
    }

    n+=1;
    count += 1;
    if(err!=4){
      gsl_histogram_increment(hN,count);
      count = 0;
      continue;
    }

    /* printing to histogram */

    gsl_histogram_increment(hG ,vG);
    gsl_histogram_increment(hRc,vR);
    gsl_histogram_increment(ha ,va);
    gsl_histogram_increment(hb ,vb);

    /* Preparing for printing */
  }while(!feof(totalVout));

  if(totalVin!=NULL) fclose(totalVin);
  if(totalVout!=NULL) fclose(totalVout);

  /**********************************************************/

  dbgPrint(11,0);

  sprintf(filename,"%s/histoInG-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInRc-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,iRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoIna-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ia,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoInb-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ib,"%f","%f");
  fclose(dadosout);

  sprintf(filename,"%s/histoOuG-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hG,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuRc-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hRc,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOua-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,ha,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOub-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hb,"%f","%f");
  fclose(dadosout);  
  sprintf(filename,"%s/histoOuN-%s.txt",folder,tag); 
  dadosout=fopen(filename,"w");
  gsl_histogram_fprintf(dadosout,hN,"%f","%f");
  fclose(dadosout);

  dbgPrint(11,0);

  sprintf(filename,"gnuplot_script.gnu");
  err=writeGnuplotScript(filename,folder,tag,nRuns,nVortex);
  if(err!=0){printf("Error printing gnuplot script\n");return err;}

  dbgPrint(12,0);

  /* histogram free - begin */
  gsl_histogram_free(hG);
  gsl_histogram_free(hRc);
  gsl_histogram_free(ha);
  gsl_histogram_free(hb);
  gsl_histogram_free(hN);

  gsl_histogram_free(iG);
  gsl_histogram_free(iRc);
  gsl_histogram_free(ia);
  gsl_histogram_free(ib);
  /* histogram free - end*/

  dbgPrint(13,0);

  return 0;
}