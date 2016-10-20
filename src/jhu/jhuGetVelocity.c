#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "turblib.h"

int main(int argc,char** argv){
  int i,j,k,n,N=10,Ntimes=1,Nx,Ny,Nz,dN;
  char *authtoken = "com.gmail.jhelsas-b854269a";
  char *dataset = "channel";
  float time0 = 0.0F,dt=0.0065,t,time1=6.5;//0.364F;
  enum SpatialInterpolation spatialInterp = Lag6;
  enum TemporalInterpolation temporalInterp = NoTInt;
  float x0=-1.0, xf=1.0, y,y0,yf,z0,zf,dx,dy,dz;
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

  Nx = 16;
  Ny = 16;
  Nz = 1;

  N = Nx * Ny * Nz;

  float (*position)[3] = malloc(N*sizeof(*position));
  float (*velocity)[3] = malloc(N*sizeof(*velocity));  
  float (*avgVelocity)[3] = malloc((N/Nz)*sizeof(*avgVelocity));  
    
  y0=-1+0.05; yf=-1+0.5;
  x0=0;       xf=M_PI/4;
  z0=0;       zf=M_PI;
  dN = 10;
  dx = (xf-x0)/((float) Nx);
  dy = (yf-y0)/((float) Ny);
  if(Nz > 1)
    dz = (zf-z0)/((float) (Nz-1));
  else
    dz = 0.;
  
  for(k=0;k<Nz;k+=1)
    for(i=0;i<Ny;i+=1)
      for(j=0;j<Nx;j+=1){
        position[k*Ny*Nx+(i*Nx+j)][0] = x0+((float) j)*dx;
        position[k*Ny*Nx+(i*Nx+j)][1] = y0+((float) i)*dy;
        position[k*Ny*Nx+(i*Nx+j)][2] = z0+((float) k)*dz;
      }

  {
    FILE *dadosAxis = fopen("Yaxis.dat","w");
    for(i=0;i<Ny;i+=1)
     fprintf(dadosAxis,"%lf\n",y0+((float) i)*dy);
    fclose(dadosAxis); 

    dadosAxis = fopen("Xaxis.dat","w");
    for(j=0;j<Nx;j+=1)
      fprintf(dadosAxis,"%lf\n",x0+((float) j)*dx);
    fclose(dadosAxis);

    dadosAxis = fopen("axis.csv","w");
    for(i=0;i<N;i+=1)
      fprintf(dadosAxis,"%f,%f,%f\n",position[Nz*0+i][0]
                                    ,position[Nz*0+i][1]
                                    ,position[Nz*0+i][2]);
    fclose(dadosAxis);
  }

  /* Initialize gSOAP */
  soapinit();

  /* Enable exit on error.  See README for details. */
  turblibSetExitOnError(1);
  
  for(i=0;i<(N/Nz);i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] = 0.;

  for(n=0;n<Ntimes;n+=dN){
    t=time0+n*dt;
    printf("\nRequesting velocity at t=%f\n",t);
    getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
                 N, position, velocity);

    for(k=0;k<Nz;k+=1)
      for(i=0;i<(N/Nz);i+=1)
        for(j=0;j<3;j+=1)
          avgVelocity[i][j] += velocity[k*Ny*Nx+i][j];

    //for(i=0;i<N;i+=1)
    for(k=0;k<Nz;k+=1){
      for(i=0;i<(N/Nz);i+=1)
        fprintf(uFile,"%f %f\n",velocity[k*Ny*Nx+i][0],velocity[k*Ny*Nx+i][1]);
      fprintf(uFile,"\n");
    }
  }

  for(i=0;i<(N/Nz);i+=1)
    for(j=0;j<3;j+=1)
      avgVelocity[i][j] /= Ntimes;
  
  for(i = 0; i < (N/Nz); i++)
    fprintf(dadosout,"%f %f %f %f %f\n",position[Nz*0+i][0],position[Nz*0+i][1],avgVelocity[i][0]
                                       ,avgVelocity[i][1],avgVelocity[i][2]);

  fclose(dadosout);
  fclose(uFile);

  /* Free gSOAP resources */
  soapdestroy();
  free(position);
  free(velocity);

  return 0;
}
