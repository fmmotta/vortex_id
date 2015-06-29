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
  const char *genFile,*tag;
} configVar;

static int vortexIdHandler(void* user, const char* section, 
                           const char* name,const char* value)
{
  configVar *vConfig = (configVar*)user;
  
    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("Grid_Parameters", "Width")) {
        vconfig->Width = atoi(value);
    } else if (MATCH("Grid_Parameters", "Heigth")) {
        vconfig->Heigth = atoi(value);
    } else if (MATCH("user", "email")) {
        vconfig->email = strdup(value);
    } else {
        return 0;  /* unknown section/name, error */
    }
    
    return 1; /* returns sucess */
}
