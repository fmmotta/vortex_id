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
  char filename[200+1];
  float time0 = 0.0F,dt=0.0065,t,time1=6.5;//0.364F;
  enum SpatialInterpolation spatialInterp = Lag6;
  enum TemporalInterpolation temporalInterp = NoTInt;
  float x0=-1.0, xf=1.0, y,y0,yf,z0,zf,dx,dy,dz,uAvg[3],gradAvg[9];
  FILE *planeFile;

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
  float (*gradient)[9] = malloc(N*sizeof(*gradient));  
  double *data = (double*) malloc(Nx*sizeof(double));
  double (*Espect)[3] = malloc(Nx*sizeof(*Espect));
  double (*Espect9)[9] = malloc(Nx*sizeof(*Espect9));
  double (*Espect18)[18] = malloc(Nx*sizeof(*Espect18));

  gsl_fft_real_wavetable * real;
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
  
  for(j=0;j<Nx;j+=1){
    int c;

    data[j]=0.;
    Espect[j][0]=0.;
    Espect[j][1]=0.;
    Espect[j][2]=0.;

    Espect9[j][0]=0.;
    Espect9[j][1]=0.;
    Espect9[j][2]=0.;
    Espect9[j][3]=0.;
    Espect9[j][4]=0.;
    Espect9[j][5]=0.;
    Espect9[j][6]=0.;
    Espect9[j][7]=0.;
    Espect9[j][8]=0.;

    for(c=0;c<18;c+=1)
      Espect18[j][c] = 0.;
  }

  i=0.;
  planeFile=NULL;
  for(n=0;n<Ntimes;n+=1){
    t=time0+n*dt;
    
    // ---------------------------------------------------------

    printf("\nRequesting velocity at t=%f\n",t);
    sprintf(filename,"plane-%d.dat",n);
    planeFile=fopen(filename,"r");
    if(planeFile==NULL){
      printf("Downloading\n");
      getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
                   N, position, velocity);
      planeFile=fopen(filename,"w");
      if(planeFile==NULL)
        return -1;
      printf("Printing\n");
      for(k=0;k<Nz;k+=1)
        for(j=0;j<Nx;j+=1)
          fprintf(planeFile,"%lf %lf %lf\n",velocity[k*Nx+j][0]
                                           ,velocity[k*Nx+j][1]
                                           ,velocity[k*Nx+j][2]);
      fclose(planeFile);
    }
    else{
      printf("Reading\n");
      for(k=0;k<Nz;k+=1)
        for(j=0;j<Nx;j+=1)
          fscanf(planeFile,"%f %f %f\n",&(velocity[k*Nx+j][0])
                                       ,&(velocity[k*Nx+j][1])
                                       ,&(velocity[k*Nx+j][2]));

      fclose(planeFile); planeFile=NULL;
    }
    printf("Velocity retrieved\n");

    // ---------------------------------------------------------

    printf("\nRequesting velocity gradient at t=%f\n",t);
    sprintf(filename,"planeGrad-%d.dat",n);
    planeFile=fopen(filename,"r");
    if(planeFile==NULL){
      printf("Downloading\n");
      //getVelocity (authtoken, dataset, t, spatialInterp, temporalInterp, 
      //             N, position, velocity);
      getVelocityGradient (authtoken, dataset, t, FD4NoInt, temporalInterp, 
                           N, position, gradient);
      planeFile=fopen(filename,"w");
      if(planeFile==NULL)
        return -1;
      printf("Printing\n");
      for(k=0;k<Nz;k+=1)
        for(j=0;j<Nx;j+=1)
          fprintf(planeFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
                                           ,gradient[k*Nx+j][0]
                                           ,gradient[k*Nx+j][1]
                                           ,gradient[k*Nx+j][2]
                                           ,gradient[k*Nx+j][3]
                                           ,gradient[k*Nx+j][4]
                                           ,gradient[k*Nx+j][5]
                                           ,gradient[k*Nx+j][6]
                                           ,gradient[k*Nx+j][7]
                                           ,gradient[k*Nx+j][8]);
      fclose(planeFile);
    }
    else{
      printf("Reading\n");
      for(k=0;k<Nz;k+=1)
        for(j=0;j<Nx;j+=1)
          fscanf(planeFile,"%f %f %f %f %f %f %f %f %f\n"
                                       ,&(gradient[k*Nx+j][0])
                                       ,&(gradient[k*Nx+j][1])
                                       ,&(gradient[k*Nx+j][2])
                                       ,&(gradient[k*Nx+j][3])
                                       ,&(gradient[k*Nx+j][4])
                                       ,&(gradient[k*Nx+j][5])
                                       ,&(gradient[k*Nx+j][6])
                                       ,&(gradient[k*Nx+j][7])
                                       ,&(gradient[k*Nx+j][8]));

      fclose(planeFile); planeFile=NULL;
    }
    printf("Velocity gradient retrieved\n");

    // ---------------------------------------------------------
    
    uAvg[0]=uAvg[1]=uAvg[2]=0.;
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Nx;j+=1){
        uAvg[0] += velocity[k*Nx+j][0];
        uAvg[1] += velocity[k*Nx+j][1];
        uAvg[2] += velocity[k*Nx+j][2];
      }
    uAvg[0] /= Nx*Nz;
    uAvg[1] /= Nx*Nz;
    uAvg[2] /= Nx*Nz;

    // ---------------------------------------------------------
    
    for(k=0;k<9;k+=1)
      gradAvg[k]=0.;
    for(k=0;k<Nz;k+=1)
      for(j=0;j<Nx;j+=1){
        gradAvg[0] += gradient[k*Nx+j][0];
        gradAvg[1] += gradient[k*Nx+j][1];
        gradAvg[2] += gradient[k*Nx+j][2];
        gradAvg[3] += gradient[k*Nx+j][3];
        gradAvg[4] += gradient[k*Nx+j][4];
        gradAvg[5] += gradient[k*Nx+j][5];
        gradAvg[6] += gradient[k*Nx+j][6];
        gradAvg[7] += gradient[k*Nx+j][7];
        gradAvg[8] += gradient[k*Nx+j][8];
      }

    for(k=0;k<9;k+=1){
      gradAvg[k] /= Nx*Nz;
      printf("grad %d = %f\n",k,gradAvg[k]);
    }
    // ---------------------------------------------------------
    
    {
      int c;
      
      for(c=0;c<3;c+=1){
        for(k=0;k<Nz;k+=1){
          for(j=0;j<Nx;j+=1)
            data[j] = (double) (velocity[k*Nx+j][c] - uAvg[c]);

          gsl_fft_real_transform(data,1,Nx,real, work);
          for(j=0;j<(Nx+1)/2;j+=1)
            Espect[j][c] += data[2*j]*data[2*j]+data[2*j+1]*data[2*j+1];
        }
      }
    }

    // ---------------------------------------------------------

    {
      int c;
      
      for(c=0;c<9;c+=1){
        for(k=0;k<Nz;k+=1){
          for(j=0;j<Nx;j+=1)
            data[j] = (double) (gradient[k*Nx+j][c] - gradAvg[c]);

          gsl_fft_real_transform(data,1,Nx,real, work);
          for(j=0;j<(Nx+1)/2;j+=1)
            Espect9[j][c] += data[2*j]*data[2*j]+data[2*j+1]*data[2*j+1];
        }
      }
    }

    // ---------------------------------------------------------
    
  }

  for(j=0;j<(Nx+1)/2;j+=1){
    int c;
    for(c=0;c<3;c+=1)
      Espect[j][c] /= Nz*Ntimes*(xf-x0);

    for(c=0;c<9;c+=1)
      Espect9[j][c] /= Nz*Ntimes*(xf-x0);
  }

  {
    int c;
    double norm=0.;
    norm = 0.;
    for(j=0;j<(Nx+1)/2;j+=1)
      norm += Espect[j][0];

    for(c=0;c<3;c+=1)
      for(j=0;j<(Nx+1)/2;j+=1)
        Espect[j][c] /= norm;

    //norm = 0.;
    //for(j=0;j<(Nx+1)/2;j+=1)
    //  norm += Espect9[j][0];

    for(c=0;c<9;c+=1)
      for(j=0;j<(Nx+1)/2;j+=1)
        Espect9[j][c] /= norm;
  }

  {
    FILE *spctOut = fopen("spectre.dat","w");
    for(j=0;j<(Nx+1)/2;j+=1)
      fprintf(spctOut,"%lf %.14lf %.14lf %.14lf\n",(2*M_PI*j)/(xf-x0),
                                Espect[j][0],Espect[j][1],Espect[j][2]);
    fclose(spctOut);
  }

  {
    FILE *spctOut = fopen("spectreGrad.dat","w");
    for(j=0;j<(Nx+1)/2;j+=1)
      fprintf(spctOut,"%lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf"
                         " %.14lf %.14lf\n",
                  (2*M_PI*j)/(xf-x0),Espect9[j][0],Espect9[j][1],Espect9[j][2]
                                    ,Espect9[j][3],Espect9[j][4],Espect9[j][5]
                                    ,Espect9[j][6],Espect9[j][7],Espect9[j][8]);
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
