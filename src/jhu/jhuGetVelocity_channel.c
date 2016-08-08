#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "turblib.h"

int main(int argc,char** argv){
  int i,j,n,N=10,Ntimes=4000,attempts=0;
  char *authtoken = "br.ufrj.if.foster-2377d412";
  char *dataset = "channel";
  float time0 = 0.0F,dt=0.0065,t,time1=6.5;//0.364F;
  enum SpatialInterpolation spatialInterp = Lag6;
  enum TemporalInterpolation temporalInterp = NoTInt;
  float x0=-1.0, xf=1.0, y, z0, zf;
  FILE *dadosout,*uFile;

	FILE *tposition = NULL;
	tposition = fopen ("teste_position.dat","w");

  if(argc!=5){
    printf("wrong number of arguments\n");
    return 1;
  }

  N = atof(argv[1]);

  uFile = fopen(argv[2],"w");
  dadosout=fopen(argv[3],"w");
  if(dadosout==NULL||uFile==NULL){
    printf("could not open file\n");
    return 2;
  }

  y = atof(argv[4]);

  // Ny = 2048
  // Nx = 25736 = 1+(int)(8 x pi x 1024)
  // Nz = 9651  = 1+(int)(3 x pi x 1024)

  float (*position)[3] = malloc(N*sizeof(*position));
  float (*velocity)[3] = malloc(N*sizeof(*velocity));  
  float (*avgVelocity)[3] = malloc(N*sizeof(*avgVelocity));  
  /*
  x0=-1.0; xf=1.0;
  for(i=0;i<N;i+=1){
  	position[i][0] = 0.;
  	position[i][1] = x0+(((float) i)/((float) N-1))*(xf-x0);
  	position[i][2] = 0.;
  }
  */

//Para puxar a linha em Z
  z0=0; zf=3*M_PI;
  for(i=0;i<N;i+=1){
    position[i][0] = 4*M_PI;
    position[i][1] = -1+y*(0.0010006);
    position[i][2] = z0+(((float) i)/((float) N-1))*(zf-z0);
		fprintf(tposition, "%f %f\n",position[i][2],position[i][1]);  
	}



// Para puxar a linha em X  
/*
  x0=0; xf=8*M_PI;
  for(i=0;i<N;i+=1){
    position[i][0] = x0+(((float) i)/((float) N-1))*(xf-x0);
    position[i][1] = -1+y*(0.0010006);
    position[i][2] = 0.;
		fprintf(tposition, "%f %f\n",position[i][0],position[i][1]);  
	}
*/
  /* Initialize gSOAP */
  soapinit();

  /* Enable exit on error.  See README for details. */
  //turblibSetExitOnError(1);
  
  for(i=0;i<N;i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] = 0.;

  for(n=0;n<Ntimes;n+=1){
    t=time0+n*dt;
    printf("\nRequesting velocity at t=%f\n",t);
    //getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp,N, position, velocity);

	  while (getVelocity (authtoken, dataset, time0, spatialInterp, temporalInterp, N, position, velocity) != SOAP_OK) {
  	  if (attempts++ > 100) {
  	    printf("Fatal Error: too many failures\n");
        exit(1);
      } 
  	  else {
    		printf("Temporary Error: %s\n", turblibGetErrorString());
		  	printf("After error %d, returning to download...\n",attempts);
			  sleep(60);      
			}
		}

    for(i=0;i<N;i+=1)
      for(j=0;j<3;j+=1)
        avgVelocity[i][j] += velocity[i][j];

    for(i=0;i<N;i+=1)
      fprintf(uFile,"%f %f %f %f\n",position[i][2],velocity[i][0]
                                   ,velocity[i][1],velocity[i][2]);
    fprintf(uFile,"\n");
  }

  for(i=0;i<N;i+=1)
    fprintf(dadosout,"%f %f %f %f\n",position[i][2],velocity[i][0]
                                   ,avgVelocity[i][1],avgVelocity[i][2]);
  for (i = 0; i < N; i++)
    fprintf(dadosout,"%f %f %f %f\n", z0+(((float) i)/((float) N-1))*(zf-z0), 
                     avgVelocity[i][0],avgVelocity[i][1],avgVelocity[i][2]);

  for(i=0;i<N;i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] /= Ntimes;


  fclose(dadosout);
  fclose(uFile);

  /* Free gSOAP resources */
  soapdestroy();
  free(position);
  free(velocity);

  return 0;
}
