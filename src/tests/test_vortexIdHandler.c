#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "inputManager.h"
#include "ini.h"

int main(int argc, char* argv[])
{
  int err;
  configVar config;
   
  err=initConfig(&config);

  if (ini_parse("cfg/test.ini", vortexIdHandler, &config) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
    
  err=printConfig(&config);

  err=freeConfig(&config);
 
  return 0;
}