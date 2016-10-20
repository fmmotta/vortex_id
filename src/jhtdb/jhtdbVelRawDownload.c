#include <math.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "ini.h"
#include "turblib.h"
#include "inputManager.h"

#define DEBUG_PRINT false

#define dbgPrint(num,num2) if(DEBUG_PRINT) printf("check point - %d-%d\n",(num),(num2))

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));      \
                                  if(ptr==NULL){                               \
                                    printf("memory not allocked\n");           \
                                    return 1;                                  \
                                  }                                            \
                                  else{                                        \
                                    int iterIndx;                              \
                                    for(iterIndx=0;iterIndx<(size);iterIndx+=1)\
                                      ptr[iterIndx]=(type) 0;                  \
                                  }                                            \

int main(int argc,char **argv){
  configVar cfg;
  FILE *dataout;
  float *rawData;
  float t0,dt,tf,t;
  int iT=0,Tw=0,n,err=0,p;
  int iX=0,iY=0,iZ=0,Xw=0,Yw=0,Zw=0;
  char filename[400+1],authtoken[400+1];
  char dataset[400+1],jhtdbFolder[400+1];
  
  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  err=initConfig(&cfg);
  
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0){
    printf("Can't load .ini file\n");
    return 1;
  }

  iT = cfg.jhtdb_iT; Tw = cfg.jhtdb_Tw;
  iX = cfg.jhtdb_Raw_iX; iY = cfg.jhtdb_Raw_iY; iZ = cfg.jhtdb_Raw_iZ;
  Xw = cfg.jhtdb_Raw_Xw; Yw = cfg.jhtdb_Raw_Yw; Zw = cfg.jhtdb_Raw_Zw;
  t0 = (float)cfg.jhtdb_t0; tf = (float)cfg.jhtdb_tf;
  strcpy(dataset,cfg.jhtdb_dataset);
  strcpy(jhtdbFolder,cfg.jhtdb_folder); 
  strcpy(authtoken,cfg.jhtdb_authToken);

  if(DEBUG_PRINT){
    err = printConfig(&cfg);
    if(err!=0){printf("problems with printConfig\n"); return err;}
  }
  
  printf("N = %d\n",Xw*Yw*Zw);
  fieldAlloc(rawData,3*Xw*Yw*Zw,float);
  
  err=mkdir(jhtdbFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %s - %d\n",jhtdbFolder,err);
    if(rawData!=NULL){free(rawData);rawData=NULL;}
    return err;
  }

  dbgPrint(3,0);
  t0 = 0.0F;
  dt = 0.0065F;
  
  dbgPrint(2,0);
  /* Initialize gSOAP */
  soapinit();

  dbgPrint(2,1);
  /* Enable exit on error.  See README for details. */
  turblibSetExitOnError(1);
  dbgPrint(2,2);
 
  for(n=iT;n<iT+Tw;n+=1){
    t = t0+dt*((float)n);
    dbgPrint(3,1);
    printf("Requesting raw velocity data: t=%lf (n=%d)\n",t,n);
    err=0;
    dbgPrint(3,2);
    if(DEBUG_PRINT){
      printf("---------------------------\n");
      printf("authtoken=%s\n",authtoken);
      printf("dataset=%s\n",dataset);
      printf("t=%lf\n",t);
      printf("i = (%d,%d,%d)\n",iX,iY,iZ);
      printf("W = (%d,%d,%d)\n",Xw,Yw,Zw);
      printf("rawData=%p\n",rawData);
      printf("---------------------------\n");
    }
    if(DEBUG_PRINT){
      printf("authtoken = %s\n",authtoken);
      printf("dataset = %s\n",dataset);
      printf("t = %lf\n",t);
      printf("start = (%d,%d,%d)\n",iX,iY,iZ);
      printf("width = (%d,%d,%d)\n",Xw,Yw,Zw);
    }
    err=getRawVelocity(authtoken,dataset,t,iX,iY,iZ,Xw,Yw,Zw,(char*)rawData);
    dbgPrint(3,3);
    if(err!=0){
      printf("Problem executing getRawVelocity: t=%lf (n=%d)\n",t,n);
      break;
    }
    dbgPrint(3,4);
    sprintf(filename,"%s/slice-%d.dat",jhtdbFolder,n);
    dbgPrint(3,5);
    dataout=fopen(filename,"w");
    dbgPrint(3,5);
    if(dataout==NULL){
      printf("Couldn't open output file - %s\n",filename);
      break;
    }
    dbgPrint(3,6);

    for(p=0;p<Xw*Yw*Zw;p+=1)
      fprintf(dataout,"%lf %lf %lf\n",rawData[3*p+0],rawData[3*p+1]
                                     ,rawData[3*p+2]);

    dbgPrint(3,7);
    fclose(dataout); dataout=NULL;
    dbgPrint(3,8);
  }

  /* Free gSOAP resources */
  soapdestroy();
 
  if(rawData!=NULL){free(rawData);rawData=NULL;}

  return 0;
}