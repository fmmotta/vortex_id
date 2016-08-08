#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "turblib.h"

int main(int argc,char** argv){
  int i,j,n,N=10,Ntimes=4000,Nx,Ny,Nz;
  char *authtoken = "com.gmail.jhelsas-b854269a";
  char *dataset = "channel";
  float time0 = 0.0F,dt=0.0065,t,time1=6.5;//0.364F;
  enum SpatialInterpolation spatialInterp = Lag6;
  enum TemporalInterpolation temporalInterp = NoTInt;
  float x0=-1.0, xf=1.0, y,y0,yf;
  FILE *dadosout,*uFile;

  if(argc!=3){
    printf("wrong number of arguments\n");
    return 1;
  }

  uFile = fopen(argv[1],"w");
  dadosout=fopen(argv[2],"w");
  if(dadosout==NULL||uFile==NULL){
    printf("could not open file\n");
    return 2;
  }

  Nx = 301;
  Ny = 96;

  N = Nx * Ny;

  float (*position)[3] = malloc(N*sizeof(*position));
  float (*velocity)[3] = malloc(N*sizeof(*velocity));  
  float (*avgVelocity)[3] = malloc(N*sizeof(*avgVelocity));  
    
  y0=-0.75; yf = -0.5;
  x0=0; xf=M_PI/4;
  for(i=0;i<Ny;i+=1)
    for(j=0;j<Nx;j+=1){
      position[i*Nx+j][0] = x0+(((float) j)/((float) (Nx-1)))*(xf-x0);
      position[i*Nx+j][1] = y0+(((float) i)/((float) (Ny-1)))*(yf-y0);
      position[i*Nx+j][2] = 0.;
    }

  {
    FILE *dadosAxis = fopen("Yaxis.dat","w");
    for(i=0;i<Ny;i+=1)
      fprintf(dadosAxis,"%lf\n",y0+(((float) i)/((float) (Ny-1)))*(yf-y0));
    fclose(dadosAxis); 

    dadosAxis = fopen("Xaxis.dat","w");
    for(j=0;j<Nx;j+=1)
      fprintf(dadosAxis,"%lf\n",x0+(((float) j)/((float) (Nx-1)))*(xf-x0));
    fclose(dadosAxis);

    dadosAxis = fopen("axis.csv","w");
    for(i=0;i<N;i+=1)
      fprintf(dadosAxis,"%f,%f,%f\n",position[i][0],position[i][1],position[i][2]);
    fclose(dadosAxis);
  }

  /* Initialize gSOAP */
  soapinit();

  /* Enable exit on error.  See README for details. */
  turblibSetExitOnError(1);
  
  for(i=0;i<N;i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] = 0.;

  for(n=0;n<Ntimes;n+=1){
    t=time0+n*dt;
    printf("\nRequesting velocity at t=%f\n",t);
    getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
                 N, position, velocity);

    for(i=0;i<N;i+=1)
      for(j=0;j<3;j+=1)
        avgVelocity[i][j] += velocity[i][j];

    for(i=0;i<N;i+=1)
      fprintf(uFile,"%f %f\n",velocity[i][0],velocity[i][1]);
    fprintf(uFile,"\n");
  }

  for(i=0;i<N;i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] /= Ntimes;
  
  for (i = 0; i < N; i++)
    fprintf(dadosout,"%f %f %f %f %f\n",position[i][0],position[i][1],avgVelocity[i][0]
                                       ,avgVelocity[i][1],avgVelocity[i][2]);

  fclose(dadosout);
  fclose(uFile);

  /* Free gSOAP resources */
  soapdestroy();
  free(position);
  free(velocity);

  return 0;
}
