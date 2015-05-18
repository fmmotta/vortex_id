#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mt64.h"

float distf(float a,float b){
  return sqrt(a*a+b*b);
}

int genLOseenUniformList(float Gmin,float Gmax, float rmin,float rmax,
	                     float *xmin,float *xmax, unsigned long long int seed,
	                     int numVortex,float **parVortexOut){
  int i,j,taint=0;
  static stat_counter=0;
  float G,rc,a,b,r,*parVortex;

  if(fabs(xmax[0]-xmin[0])<=0)
    return 1;
  if(fabs(xmax[1]-xmin[1])<=0)
    return 2;

  if(stat_counter==0){
    stat_counter=1;
    init_genrand64(seed);
  }
  
  if(*parVortexOut==NULL)
    parVortex=(float*)malloc(4*numVortex*sizeof(float));
  else
    parVortex=*parVortexOut;

  for(i=0;i<numVortex;i+=1){

    do{
      G  = Gmin+(Gmax-Gmin)*genrand64_real1();
      rc = rmin+(rmax-rmin)*genrand64_real1();
      a  = xmin[0]+(xmax[0]-xmin[0])*genrand64_real1();
      b  = xmin[1]+(xmax[1]-xmin[1])*genrand64_real1();
     
      taint=0;
      for(j=0;j<i;j+=1){
        r = distf(a-parVortex[4*j+2],b-parVortex[4*j+3]);
        if(r< 1.2*(rc+parVortex[4*j+1])) /* Magic Number: change later*/
          taint=1;
      }
    }while(taint!=0);

    parVortex[4*i+0]=G;
    parVortex[4*i+1]=rc;
    parVortex[4*i+2]=a;
    parVortex[4*i+3]=b;
  }

  *parVortexOut=parVortex;

  return 0;
}

int genLOseenBinaryList(float Gmin,float Gmax, float rmin,float rmax,
	                    float *xmin,float *xmax, unsigned long long int seed,
	                    int numVortex,float **parVortexOut){
  int i,j,taint=0;
  static stat_counter=0;
  float G,rc,a,b,r,*parVortex;

  if(fabs(xmax[0]-xmin[0])<=0)
    return 1;
  if(fabs(xmax[1]-xmin[1])<=0)
    return 2;

  if(stat_counter==0){
    stat_counter=1;
    init_genrand64(seed);
  }

  if(*parVortexOut==NULL)
    parVortex=(float*)malloc(4*numVortex*sizeof(float));
  else
    parVortex=*parVortexOut;
  
  for(i=0;i<numVortex;i+=1){
    do{
      G  = Gmin+(Gmax-Gmin)*((float)(genrand64_int64()%2));
      rc = rmin+(rmax-rmin)*((float)(genrand64_int64()%2));
      a  = xmin[0]+(xmax[0]-xmin[0])*(genrand64_real1());
      b  = xmin[1]+(xmax[1]-xmin[1])*(genrand64_real1());
     
      taint=0;
      for(j=0;j<i;j+=1){
        r = distf(a-parVortex[4*j+2],b-parVortex[4*j+3]);
        if(r< 1.2*(rc+parVortex[4*j+1])) /* Magic Number: change later*/
          taint=1;
      }
    }while(taint!=0);

    parVortex[4*i+0]=G;
    parVortex[4*i+1]=rc;
    parVortex[4*i+2]=a;
    parVortex[4*i+3]=b;
  }

  *parVortexOut=parVortex;

  return 0;
}

