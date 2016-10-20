
typedef struct configVar {
  long long int seed;
  int Width,Height,nVortex,nRuns;
  int runType,genType,dim,adaptive,calcMode;
  int numG,numRc,hNG,hNRc,hNa,hNb,hNN,hNmax;
  int Nx,Ny,Nz,pType,pIndex,Nsnapshots,planeNum,*pln;
  double Gmin,Gmax,RcMin,RcMax;
  double *Glist,*Rclist,x0[2],xf[2];
  double swThresh, sndSwThresh,cutoff,v0y0;
  double hGmin,hGmax,hRcMin,hRcMax,t0,dt;
  char *genFile,*tag,*folder,*FOAMfolder,*bkgFile;
  
  int jhtdb_iT,jhtdb_Tw;
  int jhtdb_Raw_iX,jhtdb_Raw_iY,jhtdb_Raw_iZ;
  int jhtdb_Raw_Xw,jhtdb_Raw_Yw,jhtdb_Raw_Zw;
  double jhtdb_t0,jhtdb_dt,jhtdb_tf;
  double jhtdb_x0,jhtdb_y0,jhtdb_z0;
  double jhtdb_dx,jhtdb_dy,jhtdb_dz;
  char *jhtdb_authToken,*jhtdb_dataset,*jhtdb_folder;
} configVar;

int initConfig(configVar *cfg);

int freeConfig(configVar *cfg);

int vortexIdHandler(void* user, const char* section, 
                    const char* name,const char* value);

int printConfig(configVar *cfg);
