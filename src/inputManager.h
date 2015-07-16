
typedef struct configVar {
  long long int seed;
  int Width,Height,nVortex,nRuns;
  int runType,genType,dim,adaptive;
  int numG,numRc,hNG,hNRc,hNa,hNb,hNN,hNmax;
  float Gmin,Gmax,RcMin,RcMax;
  float *Glist,*Rclist,x0[2],xf[2];
  float swThresh, sndSwThresh,cutoff,v0y0;
  float hGmin,hGmax,hRcMin,hRcMax;
  char *genFile,*tag,*folder;
} configVar;

int initConfig(configVar *cfg);

int freeConfig(configVar *cfg);

int vortexIdHandler(void* user, const char* section, 
                    const char* name,const char* value);

int printConfig(configVar *cfg);