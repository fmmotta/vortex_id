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
  float t0,dt,tf,t;
  float x0,y0,z0,xf,yf,zf,dx,dy,dz;
  int iT=0,Tw=0,n,err=0,p,N,i,j,k;
  int iX=0,iY=0,iZ=0,Xw=0,Yw=0,Zw=0;
  char filename[400+1],authtoken[400+1];
  char dataset[400+1],jhtdbFolder[400+1];
  enum SpatialInterpolation spatialInterp = NoSInt;
  enum TemporalInterpolation temporalInterp = NoTInt;
  
  dbgPrint(0,0);
  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  dbgPrint(1,0);
  err=initConfig(&cfg);
  
  dbgPrint(2,0);
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0){
    printf("Can't load .ini file\n");
    return 1;
  }

  dbgPrint(3,0);
  iT = cfg.jhtdb_iT; Tw = cfg.jhtdb_Tw;
  iX = cfg.jhtdb_Raw_iX; iY = cfg.jhtdb_Raw_iY; iZ = cfg.jhtdb_Raw_iZ;
  Xw = cfg.jhtdb_Raw_Xw; Yw = cfg.jhtdb_Raw_Yw; Zw = cfg.jhtdb_Raw_Zw;
  t0 = (float)cfg.jhtdb_t0; tf = (float)cfg.jhtdb_tf;
  x0 = (float)cfg.jhtdb_x0; dx = (float)cfg.jhtdb_dx;
  y0 = (float)cfg.jhtdb_y0; dy = (float)cfg.jhtdb_dy;
  z0 = (float)cfg.jhtdb_z0; dz = (float)cfg.jhtdb_dz;
  strcpy(dataset,cfg.jhtdb_dataset);
  strcpy(jhtdbFolder,cfg.jhtdb_folder); 
  strcpy(authtoken,cfg.jhtdb_authToken);

  printf("%f %f %f\n",dx,dy,dz);

  dbgPrint(4,0);
  N=Xw*Yw*Zw;
  float (*position)[3] = malloc(N*sizeof(*position));
  float (*velocity)[3] = malloc(N*sizeof(*velocity));

  dbgPrint(5,0);
  //dx = (xf-x0)/((float) Xw);
  //dy = (yf-y0)/((float) Yw);
  //dz = 0.0;
  xf = x0+Xw*dx;
  yf = y0+Yw*dy;
  zf = z0+Zw*dz;

  dbgPrint(6,0);
  for(k=0;k<Zw;k+=1)
    for(i=0;i<Yw;i+=1)
      for(j=0;j<Xw;j+=1){
        position[k*Yw*Xw+(i*Xw+j)][0] = x0+((float) j)*dx;
        position[k*Yw*Xw+(i*Xw+j)][1] = y0+((float) i)*dy;
        position[k*Yw*Xw+(i*Xw+j)][2] = z0+((float) k)*dz;
      }

  dbgPrint(7,0);
  {
    FILE *dadosAxis;
    sprintf(filename,"%s/Yaxis.dat",jhtdbFolder);
    dadosAxis = fopen(filename,"w");
    for(p=0;p<Yw;p+=1)
     fprintf(dadosAxis,"%lf\n",y0+((float) p)*dy);
    fclose(dadosAxis); 

    sprintf(filename,"%s/Xaxis.dat",jhtdbFolder);
    dadosAxis = fopen(filename,"w");
    for(p=0;p<Xw;p+=1)
      fprintf(dadosAxis,"%lf\n",x0+((float) p)*dx);
    fclose(dadosAxis); 

    sprintf(filename,"%s/Zaxis.dat",jhtdbFolder);
    dadosAxis = fopen(filename,"w");
    for(p=0;p<Zw;p+=1)
      fprintf(dadosAxis,"%lf\n",z0+((float) p)*dz);
    fclose(dadosAxis);

    sprintf(filename,"%s/axis.dat",jhtdbFolder);
    dadosAxis = fopen(filename,"w");
    for(p=0;p<N;p+=1)
      fprintf(dadosAxis,"%f,%f,%f\n",position[Zw*0+i][0]
                                    ,position[Zw*0+i][1]
                                    ,position[Zw*0+i][2]);
    fclose(dadosAxis);
  }
  dbgPrint(8,0);

  if(DEBUG_PRINT){
    err = printConfig(&cfg);
    if(err!=0){printf("problems with printConfig\n"); return err;}
  }
  
  dbgPrint(9,0);
  err=mkdir(jhtdbFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %s - %d\n",jhtdbFolder,err);
    if(position!=NULL) free(position);
    if(velocity!=NULL) free(velocity);
    return err;
  }
  dbgPrint(10,0);

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
 
  dbgPrint(11,0);
  for(n=iT;n<iT+Tw;n+=1){
    t = t0+dt*((float)n);
    dbgPrint(3,1);
    printf("Requesting interpolated velocity data: t=%lf (n=%d)\n",t,n);
    err=0;
    dbgPrint(3,2);
    if(DEBUG_PRINT){
      printf("---------------------------\n");
      printf("authtoken=%s\n",authtoken);
      printf("dataset=%s\n",dataset);
      printf("t=%lf\n",t);
      printf("i = (%d,%d,%d)\n",iX,iY,iZ);
      printf("W = (%d,%d,%d)\n",Xw,Yw,Zw);
      printf("---------------------------\n");
    }
    //err=getRawVelocity(authtoken,dataset,t,iX,iY,iZ,Xw,Yw,Zw,(char*)rawData);
    err=getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
                     N, position, velocity);
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
      fprintf(dataout,"%lf %lf %lf\n",velocity[p][0],velocity[p][1],
                                      velocity[p][2]);

    dbgPrint(3,7);
    fclose(dataout); dataout=NULL;
    dbgPrint(3,8);
  }
  dbgPrint(12,0);

  /* Free gSOAP resources */
  soapdestroy();
  dbgPrint(13,0);
 
  if(position!=NULL) free(position);
  if(velocity!=NULL) free(velocity);

  dbgPrint(14,0);
  return 0;
}