#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inih/ini.h"
#include "inputManager.h"

typedef struct {
  int Width,Height,nVortex,nRuns;
  int runType,seed,dim,adaptive;
  int numG,numRc;
  float Gmin,Gmax,RcMin,RcMax;
  float *Glist,*Rclist;
  float swThresh, sndSwThresh;
  char *genFile,*tag;
} configVar;

int initConfig(configVar *cfg){
  cfg->Width=0; cfg->Heigth=0; cfg->nVortex=0; cfg->nRuns=0;
  cfg->runType=0; cfg->seed=0; cfg->dim=2; cfg->adaptive=0;
  cfg->numG =0; cfg->numRc =0;
  cfg->Gmin=0.; cfg->Gmax =0.; cfg->RcMin=0.; cfg->RcMax=0.;
  cfg->Glist=NULL; cfg->Rclist=NULL; 
  cfg->swThresh=0.; cfg->sndSwThresh=0.;
  cfg->genFile=NULL; cfg->tag=NULL;

  return 0;
}

int freeConfig(configVar *cfg){
  cfg->Width=0; cfg->Heigth=0; cfg->nVortex=0; cfg->nRuns=0;
  cfg->runType=0; cfg->seed=0; cfg->dim=2; cfg->adaptive=0;
  cfg->numG =0; cfg->numRc =0;
  cfg->Gmin=0.; cfg->Gmax =0.; cfg->RcMin=0.; cfg->RcMax=0.;
  cfg->swThresh=0.; cfg->sndSwThresh=0.;

  if(cfg->Glist!=NULL)
    free(cfg->Glist);
  if(cfg->Rclist!=NULL)
    free(cfg->Rclist);
  if(cfg->genFile!=NULL)
    free(cfg->genFile);
  if(cfg->tag)
    free(cfg->tag);

  cfg->Glist=NULL; cfg->Rclist=NULL; 
  cfg->genFile=NULL; cfg->tag=NULL;

  return 0;  
}

int vortexIdHandler(void* user, const char* section, 
                    const char* name,const char* value)
{
  configVar *vConfig = (configVar*)user;
  
  #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
  if (MATCH("Grid_Parameters", "Width"))
    vConfig->Width = atoi(value);
  else if (MATCH("Grid_Parameters", "Heigth"))
    vConfig->Heigth = atoi(value);
  else if (MATCH("RunInfo", "Tag"))
    vConfig->tag = strdup(value);
  else if (MATCH("",""))
  else
    return 0;  /* unknown section/name, error */
    
  return 1; /* returns sucess */
}
