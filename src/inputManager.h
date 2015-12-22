
typedef struct configVar {
  long long int seed;
  int Width,Height,nVortex,nRuns;
  int runType,genType,dim,adaptive;
  int numG,numRc,hNG,hNRc,hNa,hNb,hNN,hNmax;
  int Nx,Ny,Nz,pType,pIndex;
  double Gmin,Gmax,RcMin,RcMax;
  double *Glist,*Rclist,x0[2],xf[2];
  double swThresh, sndSwThresh,cutoff,v0y0;
  double hGmin,hGmax,hRcMin,hRcMax;
  char *genFile,*tag,*folder,*FOAMfolder;
} configVar;

int initConfig(configVar *cfg);

int freeConfig(configVar *cfg);

int vortexIdHandler(void* user, const char* section, 
                    const char* name,const char* value);

int printConfig(configVar *cfg);