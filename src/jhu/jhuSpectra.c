#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "turblib.h"

int main(int argc,char** argv){
  int i,j,k,n,N=10,Ntimes=4000,Nx,Ny,Nz;
  char *authtoken = "com.gmail.jhelsas-b854269a";
  char *dataset = "channel";
  float time0 = 0.0F,dt=0.0065,t,time1=6.5;//0.364F;
  enum SpatialInterpolation spatialInterp = Lag6;
  enum TemporalInterpolation temporalInterp = NoTInt;
  float x0=-1.0, xf=1.0, y,y0,yf,z0,zf,dx,dy,dz;

  Nx = 2048;
  Ny = 1;
  Nz = 1536;
  Ntimes = 1;

  y0=0.0; yf=0.0;
  x0=0;      xf=8*M_PI;
  z0=0;      zf=3*M_PI;
  dx = (xf-x0)/((float) Nx);
  dy = (yf-y0)/((float) Ny);
  dz = (zf-z0)/((float) Nz);

  N = Nx * Ny * Nz;

  float (*position)[3] = malloc(N*sizeof(*position));
  float (*velocity)[3] = malloc(N*sizeof(*velocity));  
  double *data = (double*) malloc(Nx*sizeof(double));
  double *Espect = (double*)malloc(Nx*sizeof(double));

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  work = gsl_fft_real_workspace_alloc (Nx);
  real = gsl_fft_real_wavetable_alloc (Nx);
  
  for(i=0;i<Ny;i+=1)
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Nx;j+=1){
        int idx = i*Nx*Nz+(k*Nx+j); // k*Ny*Nx+(i*Nx+j); 
        position[idx][0] = x0+((float) j)*dx;
        position[idx][1] = y0+((float) i)*dy;
        position[idx][2] = z0+((float) k)*dz;
      }

  {
    FILE *dadosAxis = fopen("Yaxis.dat","w");
    for(i=0;i<Ny;i+=1)
     fprintf(dadosAxis,"%lf\n",y0+((float) i)*dy);
    fclose(dadosAxis); 

    dadosAxis = fopen("Zaxis.dat","w");
    for(k=0;k<Nz;k+=1)
      fprintf(dadosAxis,"%lf\n",z0+((float) k)*dz);
    fclose(dadosAxis);

    dadosAxis = fopen("Xaxis.dat","w");
    for(j=0;j<Nx;j+=1)
      fprintf(dadosAxis,"%lf\n",x0+((float) j)*dx);
    fclose(dadosAxis);

    dadosAxis = fopen("axis.csv","w");
    for(i=0;i<N;i+=1)
      fprintf(dadosAxis,"%f,%f,%f\n",position[i][0]
                                    ,position[i][1]
                                    ,position[i][2]);
    fclose(dadosAxis);
  }

  /* Initialize gSOAP */
  soapinit();

  /* Enable exit on error.  See README for details. */
  turblibSetExitOnError(1);
  
  for(i=0;i<Nx;i+=1){
    data[i]=0.;
    Espect[i]=0.;
  }

  for(n=0;n<Ntimes;n+=1){
    t=time0+n*dt;
    printf("\nRequesting velocity at t=%f\n",t);
    getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
                 N, position, velocity);
    printf("Velocity retrieved\n");

    for(k=0;k<Nz;k+=1){
      printf("Processing spectre contribution for line k=%d\n",k);
      for(j=0;j<Nx;j+=1)
        data[i] = (double) velocity[k*Nx+j][0];

      gsl_fft_real_transform (data,1,Nx,real, work);
      for(j=0;j<Nx;j+=1)
        Espect[j] += data[j]*data[j];
    }
  }

  for(j=0;j<Nx;j+=1)
    Espect[j] /= Nz*Ntimes*(xf-x0);

  {
    FILE *spctOut = fopen("spectre.dat","w");
    for(j=0;j<Nx;j+=1)
      fprintf(spctOut,"%lf %lf\n",position[j][0],Espect[j]);
    fclose(spctOut);
  }

  /* Free gSOAP resources */

  soapdestroy();
  free(position);
  free(velocity);
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);
  free(data);

  return 0;
}
