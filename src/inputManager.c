#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputManager.h"

int initConfig(configVar *cfg){
  cfg->Width=0; cfg->Height=0; cfg->nVortex=0; cfg->nRuns=0;
  cfg->runType=-1; cfg->seed=0; cfg->dim=2; cfg->adaptive=0;
  cfg->numG =0; cfg->numRc =0; cfg->hNmax=0; cfg->genType=0;
  cfg->Gmin=0.; cfg->Gmax =0.; cfg->RcMin=0.; cfg->RcMax=0.;
  cfg->hNG=0; cfg->hNRc=0; cfg->hNa=0; cfg->hNb=0; cfg->hNN=0;
  cfg->Nx=0;cfg->Ny=0;cfg->Nz=0;cfg->pType=0;cfg->pIndex=0;
  cfg->jhtdb_iT=0,cfg->jhtdb_Tw=0;
  cfg->jhtdb_Raw_iX=0,cfg->jhtdb_Raw_iY=0,cfg->jhtdb_Raw_iZ=0;
  cfg->jhtdb_Raw_Xw=0,cfg->jhtdb_Raw_Yw=0,cfg->jhtdb_Raw_Zw=0;
  cfg->jhtdb_t0=0.0,  cfg->jhtdb_dt=0.0,  cfg->jhtdb_tf=0.0;
  cfg->jhtdb_x0=0.0,  cfg->jhtdb_y0=0.0,  cfg->jhtdb_z0=0.0;
  cfg->jhtdb_dx=0.0,  cfg->jhtdb_dy=0.0,  cfg->jhtdb_dz=0.0; 
  cfg->planeNum=0; 
  cfg->calcMode=-1; cfg->FOAMfolder=NULL;
  cfg->Glist=NULL; cfg->Rclist=NULL; 
  cfg->swThresh=0.; cfg->sndSwThresh=0.; cfg->cutoff=0.;
  cfg->genFile=NULL; cfg->tag=NULL; cfg->pln=NULL;
  cfg->folder=NULL; cfg->bkgFile=NULL; 
  cfg->jhtdb_authToken=NULL,cfg->jhtdb_dataset=NULL; cfg->jhtdb_folder=NULL;

  return 0;
}

int freeConfig(configVar *cfg){
  cfg->Width=0; cfg->Height=0; cfg->nVortex=0; cfg->nRuns=0;
  cfg->runType=-1; cfg->seed=0; cfg->dim=2; cfg->adaptive=0;
  cfg->numG =0; cfg->numRc =0; cfg->hNmax=0; cfg->adaptive=0;
  cfg->Gmin=0.; cfg->Gmax =0.; cfg->RcMin=0.; cfg->RcMax=0.;
  cfg->swThresh=0.; cfg->sndSwThresh=0.; cfg->cutoff=0.;cfg->planeNum=0; 
  cfg->hNG=0; cfg->hNRc=0; cfg->hNa=0; cfg->hNb=0; cfg->hNN=0;
  cfg->Nx=0;cfg->Ny=0;cfg->Nz=0;cfg->pType=0;cfg->pIndex=0;
  cfg->Nsnapshots=0; cfg->t0 = 0.; cfg->dt=0.; cfg->calcMode-=-1;
  cfg->jhtdb_iT=0,    cfg->jhtdb_Tw=0;
  cfg->jhtdb_Raw_iX=0,cfg->jhtdb_Raw_iY=0,cfg->jhtdb_Raw_iZ=0;
  cfg->jhtdb_Raw_Xw=0,cfg->jhtdb_Raw_Yw=0,cfg->jhtdb_Raw_Zw=0;
  cfg->jhtdb_t0=0.0,  cfg->jhtdb_dt=0.0,  cfg->jhtdb_tf=0.0;
  cfg->jhtdb_x0=0.0,  cfg->jhtdb_y0=0.0,  cfg->jhtdb_z0=0.0;
  cfg->jhtdb_dx=0.0,  cfg->jhtdb_dy=0.0,  cfg->jhtdb_dz=0.0; 
  

  if(cfg->Glist!=NULL)
    free(cfg->Glist);
  if(cfg->Rclist!=NULL)
    free(cfg->Rclist);
  if(cfg->genFile!=NULL)
    free(cfg->genFile);
  if(cfg->tag!=NULL)
    free(cfg->tag);
  if(cfg->folder!=NULL)
    free(cfg->folder);
  if(cfg->pln!=0)
    free(cfg->pln);
  if(cfg->FOAMfolder!=NULL)
    free(cfg->FOAMfolder);
  if(cfg->bkgFile!=NULL)
    free(cfg->bkgFile);
  if(cfg->jhtdb_authToken!=NULL)
    free(cfg->jhtdb_authToken);
  if(cfg->jhtdb_dataset!=NULL)
    free(cfg->jhtdb_dataset);
  if(cfg->jhtdb_folder!=NULL)
    free(cfg->jhtdb_folder);

  cfg->Glist=NULL; cfg->Rclist=NULL; 
  cfg->genFile=NULL; cfg->tag=NULL;
  cfg->folder=NULL;cfg->FOAMfolder=NULL;
  cfg->pln=NULL; cfg->bkgFile=NULL;
  cfg->jhtdb_authToken=NULL;
  cfg->jhtdb_dataset=NULL;
  cfg->jhtdb_folder=NULL;

  return 0;  
}

int vortexIdHandler(void* user, const char* section, 
                    const char* name,const char* value)
{
  int i;
  char *p;
  configVar *vConfig = (configVar*)user;
  
  #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
  if (MATCH("Grid-Parameters", "Width"))
    vConfig->Width = atoi(value);
  else if (MATCH("Grid-Parameters", "Height"))
    vConfig->Height = atoi(value);
  else if (MATCH("Grid-Parameters","x0"))
    sscanf(value,"%lf %lf",&(vConfig->x0[0]), &(vConfig->x0[1]) );
  else if (MATCH("Grid-Parameters","xf"))
    sscanf(value,"%lf %lf",&(vConfig->xf[0]), &(vConfig->xf[1]) );
  else if (MATCH("Vortex-Generation","type"))
    vConfig->genType = atoi(value);
  else if (MATCH("Vortex-Generation","number-of-vortices"))
    vConfig->nVortex = atoi(value);
  else if (MATCH("Vortex-Generation","seed"))
    vConfig->seed = atoll(value);
  else if (MATCH("Vortex-Generation","numG"))
    vConfig->numG = atoi(value);
  else if (MATCH("Vortex-Generation","numRc"))
    vConfig->numRc = atoi(value);
  else if (MATCH("Vortex-Generation","Gmin"))
    vConfig->Gmin = atof(value);
  else if (MATCH("Vortex-Generation","Gmax"))
    vConfig->Gmax = atof(value);
  else if (MATCH("Vortex-Generation","RcMin"))
    vConfig->RcMin = atof(value);
  else if (MATCH("Vortex-Generation","RcMax"))
    vConfig->RcMax = atof(value);
  else if (MATCH("Vortex-Generation","Vertical-Shear"))
    vConfig->v0y0 = atof(value);
  else if (MATCH("Runtime-Info","type"))
    vConfig->runType = atoi(value);
  else if (MATCH("Runtime-Info", "Tag"))
    vConfig->tag = strdup(value);
  else if (MATCH("Runtime-Info", "Parameter-Savefile"))
    vConfig->genFile = strdup(value);
  else if (MATCH("Runtime-Info", "Runs"))
    vConfig->nRuns = atoi(value);
  else if (MATCH("Runtime-Info", "Mode"))
    vConfig->calcMode = atoi(value);
  else if (MATCH("Runtime-Info", "Dimension"))
    vConfig->dim = atoi(value);
  else if (MATCH("Runtime-Info", "Folder"))
    vConfig->folder = strdup(value);
  else if (MATCH("Runtime-Info", "Background-File"))
    vConfig->bkgFile = strdup(value);
  else if (MATCH("Reconstruction-Info", "Adaptive"))
    vConfig->adaptive = atoi(value);
  else if (MATCH("Reconstruction-Info", "Circulation-Cutoff"))
    vConfig->cutoff = atof(value);
  else if (MATCH("Reconstruction-Info", "Swirling-Strength-Threshold"))
    vConfig->swThresh = atof(value);
  else if (MATCH("Reconstruction-Info", "Second-Swirling-Strength-Threshold"))
    vConfig->sndSwThresh = atof(value);
  else if(MATCH("Histogram-Parameters","bin-G"))
    vConfig->hNG = atoi(value);
  else if(MATCH("Histogram-Parameters","bin-Rc"))
    vConfig->hNRc = atoi(value);
  else if(MATCH("Histogram-Parameters","bin-a"))
    vConfig->hNa = atoi(value);
  else if(MATCH("Histogram-Parameters","bin-b"))
    vConfig->hNb = atoi(value);
  else if(MATCH("Histogram-Parameters","bin-N"))
    vConfig->hNN = atoi(value);
  else if(MATCH("Histogram-Parameters","histogram-Gmin"))
    vConfig->hGmin = atof(value);
  else if(MATCH("Histogram-Parameters","histogram-Gmax"))
    vConfig->hGmax = atof(value);
  else if(MATCH("Histogram-Parameters","histogram-RcMin"))
    vConfig->hRcMin = atof(value);
  else if(MATCH("Histogram-Parameters","histogram-RcMax"))
    vConfig->hRcMax = atof(value);
  else if(MATCH("Histogram-Parameters","histogram-Nmax"))
    vConfig->hNmax = atoi(value);
  else if(MATCH("Vortex-Generation","Glist")){
    if(vConfig->Glist==NULL)
      vConfig->Glist=(double*)malloc((vConfig->numG)*sizeof(double));
    
    p=strtok(value," ");
    i=0;
    while(p!=NULL){
      vConfig->Glist[i] = atof(p);
      i=i+1;
      p = strtok(NULL," ");
    }
  }
  else if(MATCH("Vortex-Generation","Rclist")){
    if(vConfig->Rclist==NULL)
        vConfig->Rclist=(double*)malloc((vConfig->numRc)*sizeof(double));
    p=strtok(value," ");
    i=0;
    while(p!=NULL){
      vConfig->Rclist[i] = atof(p);
      i=i+1;
      p = strtok(NULL," ");
    }
  }
  else if(MATCH("openFOAM","folder"))
    vConfig->FOAMfolder = strdup(value);
  else if(MATCH("openFOAM","Nx"))
    vConfig->Nx=atoi(value);
  else if(MATCH("openFOAM","Ny"))
    vConfig->Ny=atoi(value);
  else if(MATCH("openFOAM","Nz"))
    vConfig->Nz=atoi(value);
  else if(MATCH("openFOAM","Plane-Type"))
    vConfig->pType=atoi(value);
  else if(MATCH("openFOAM","Plane-Index"))
    vConfig->pIndex=atoi(value);
  else if(MATCH("openFOAM","Plane-Number"))
    vConfig->planeNum=atoi(value);
  else if(MATCH("openFOAM","Plane-List")){
    if(vConfig->pln==NULL)
      vConfig->pln=(int*)malloc((vConfig->planeNum)*sizeof(int));
    
    p=strtok(value," ");
    i=0;
    while(p!=NULL){
      vConfig->pln[i] = atof(p);
      i=i+1;
      p = strtok(NULL," ");
    }
  }
  else if(MATCH("openFOAM","Nsnapshots"))
    vConfig->Nsnapshots=atoi(value);
  else if(MATCH("openFOAM","t0"))
    vConfig->t0=atof(value);
  else if(MATCH("openFOAM","dt"))
    vConfig->dt=atof(value);
  else if(MATCH("JHTDB","iT"))
    vConfig->jhtdb_iT=atoi(value);
  else if(MATCH("JHTDB","Tw"))
    vConfig->jhtdb_Tw=atoi(value);
  else if(MATCH("JHTDB","Raw_iX"))
    vConfig->jhtdb_Raw_iX=atoi(value);
  else if(MATCH("JHTDB","Raw_iY"))
    vConfig->jhtdb_Raw_iY=atoi(value);
  else if(MATCH("JHTDB","Raw_iZ"))
    vConfig->jhtdb_Raw_iZ=atoi(value);
  else if(MATCH("JHTDB","Raw_Xw"))
    vConfig->jhtdb_Raw_Xw=atoi(value);
  else if(MATCH("JHTDB","Raw_Yw"))
    vConfig->jhtdb_Raw_Yw=atoi(value);
  else if(MATCH("JHTDB","Raw_Zw"))
    vConfig->jhtdb_Raw_Zw=atoi(value);
  else if(MATCH("JHTDB","t0"))
    vConfig->jhtdb_t0=atof(value);
  else if(MATCH("JHTDB","dt"))
    vConfig->jhtdb_dt=atof(value);
  else if(MATCH("JHTDB","tf"))
    vConfig->jhtdb_tf=atof(value);
  else if(MATCH("JHTDB","x0"))
    vConfig->jhtdb_x0=atof(value);
  else if(MATCH("JHTDB","y0"))
    vConfig->jhtdb_y0=atof(value);
  else if(MATCH("JHTDB","z0"))
    vConfig->jhtdb_z0=atof(value);
  else if(MATCH("JHTDB","dx"))
    vConfig->jhtdb_dx=atof(value);
  else if(MATCH("JHTDB","dy"))
    vConfig->jhtdb_dy=atof(value);
  else if(MATCH("JHTDB","dz"))
    vConfig->jhtdb_dz=atof(value);
  else if(MATCH("JHTDB","authentication-token"))
    vConfig->jhtdb_authToken=strdup(value);
  else if(MATCH("JHTDB","dataset"))
    vConfig->jhtdb_dataset=strdup(value);
  else if(MATCH("JHTDB","folder"))
    vConfig->jhtdb_folder=strdup(value);
  else 
    return 0;  /* unknown section/name, error */

  return 1; /* returns sucess */
}


int printConfig(configVar *cfg){
  int k;
  printf("Grid Parameters:\n");
  printf("Size: %d %d\n",cfg->Height,cfg->Width);
  printf("Start: %f %f\n",cfg->x0[0],cfg->x0[1]);
  printf("End: %f %f\n",cfg->xf[0],cfg->xf[1]);

  printf("\nVortex Generation Parameters\n");
  printf("Gen Type: %d\n",cfg->runType);
  printf("Number of Vortices/Run: %d\n",cfg->nVortex);
  printf("PNRG seed: %lld\n",cfg->seed);
  if(cfg->genType==0){
    printf("G interval: [%f,%f]\n",cfg->Gmin, cfg->Gmax);
    printf("Rc interval: [%f,%f]\n",cfg->RcMin, cfg->RcMax);
  }
  if(cfg->runType==1){
    printf("Number of G samples: %d\n",cfg->numG);  
    printf("Number of Rc samples: %d\n",cfg->numRc);
    printf("G samples: ");
    for(k=0;k<cfg->numG;k+=1)
      printf("%f ",cfg->Glist[k]);
    printf("\n");  
    printf("Rc samples: ");
    for(k=0;k<cfg->numRc;k+=1)
      printf("%f ",cfg->Rclist[k]);
    printf("\n");
  }
  printf("Vertical Shear:%f\n",cfg->v0y0);

  printf("\nRuntime Info:\n");
  printf("type: %d\n",cfg->runType);
  printf("mode: %d\n",cfg->calcMode);
  printf("tag: %s\n",cfg->tag);
  printf("genFile: %s\n",cfg->genFile);
  printf("number of Runs: %d\n",cfg->nRuns);
  printf("dimension : %d\n",cfg->dim);
  printf("Results Folder : %s\n",cfg->folder);
  printf("Background File : %s\n",cfg->bkgFile);

  printf("\nReconstruction-Info:\n");
  printf("Adaptive: %d\n",cfg->adaptive);
  printf("Circulation Minimum Cutoff: %f\n",cfg->cutoff);
  printf("Swirling Strength Threshold: %f\n",cfg->swThresh);
  printf("2nd Swirling Strength Threshold: %f\n",cfg->swThresh);

  printf("\nHistogram Parameters:\n");
  printf("Number of G bins: %d\n",cfg->hNG);
  printf("Number of Rc bins: %d\n",cfg->hNRc);
  printf("Number of a bins: %d\n",cfg->hNa);
  printf("Number of b bins: %d\n",cfg->hNb);
  printf("Number of N bins: %d\n",cfg->hNN);
  printf("G Interval: [%.f , %f)\n",cfg->hGmin,cfg->hGmax);
  printf("Rc Interval: [%.f , %f)\n",cfg->hRcMin,cfg->hRcMax);
  printf("N max value: %d\n",cfg->hNmax);

  printf("\nOpenFOAM parameters:\n");
  printf("OpenFOAM folder: %s\n",cfg->FOAMfolder);
  printf("OpenFOAM Nx : %d\n",cfg->Nx);
  printf("OpenFOAM Ny : %d\n",cfg->Ny);
  printf("OpenFOAM Nz : %d\n",cfg->Nz);
  printf("Plane-Type  : %d\n",cfg->pType);
  printf("Plane-Index : %d\n",cfg->pIndex);
  printf("Plane-Number: %d\n",cfg->planeNum);
  printf("Plane-List: ");
  for(k=0;k<cfg->planeNum;k+=1)
    printf(" %d,",cfg->pln[k]);
  printf("\n");
  printf("Number of Snapshots: %d\n",cfg->Nsnapshots);
  printf("Initial Time: %f\n",cfg->t0);
  printf("Time step: %f\n",cfg->dt);

  printf("\nJHU Turbulence Database parameters:\n");
  printf("JHUTDB folder: %s\n",cfg->jhtdb_folder);
  printf("JHUTDB dataset: %s\n",cfg->jhtdb_dataset);
  printf("JHUTDB authentication-token: %s\n",cfg->jhtdb_authToken);
  printf("JHUTDB t0: %lf\n",cfg->jhtdb_t0);
  printf("JHUTDB dt: %lf\n",cfg->jhtdb_dt);
  printf("JHUTDB tf: %lf\n",cfg->jhtdb_tf);
  printf("JHUTDB iT: %d\n",cfg->jhtdb_iT);
  printf("JHUTDB Tw: %d\n",cfg->jhtdb_Tw);
  printf("JHUTDB Raw_iX: %d\n",cfg->jhtdb_Raw_iX);
  printf("JHUTDB Raw_iY: %d\n",cfg->jhtdb_Raw_iY);
  printf("JHUTDB Raw_iZ: %d\n",cfg->jhtdb_Raw_iZ);
  printf("JHUTDB Raw_Xw: %d\n",cfg->jhtdb_Raw_Xw);
  printf("JHUTDB Raw_Yw: %d\n",cfg->jhtdb_Raw_Yw);
  printf("JHUTDB Raw_Zw: %d\n",cfg->jhtdb_Raw_Zw);
  printf("JHUTDB x0: %lf\n",cfg->jhtdb_x0);
  printf("JHUTDB xf: %lf\n",cfg->jhtdb_dx);
  printf("JHUTDB y0: %lf\n",cfg->jhtdb_y0);
  printf("JHUTDB yf: %lf\n",cfg->jhtdb_dy);
  printf("JHUTDB z0: %lf\n",cfg->jhtdb_z0);
  printf("JHUTDB zf: %lf\n",cfg->jhtdb_dz);

  return 0;
}
