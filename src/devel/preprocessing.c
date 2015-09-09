#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define N 5

typedef struct openFoamIcoData {
    float ux,uy,uz,p;
} openFoamIcoData;

int main(int argc,char** argv){
  int i,j,k,l,Npre,Nu,Np;
  char buffer[1024];
  FILE *uFile,*pFile,*ouFile;
  openFoamIcoData v[N];

  if(argc!=3){
    printf("wrong number of arguments, need exactly 2 input files, velocity and pressure files");
    return 1;
  }

  Npre=20; // preamble size

  uFile = fopen(argv[1],"r");
  pFile = fopen(argv[2],"r");

  if(uFile==NULL || pFile == NULL){printf("problems opening the files\n"); return 1;}

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,uFile);
  fscanf(uFile,"%d",&Nu);

  for(i=0;i<Npre;i++)
    fgets(buffer,1024,pFile);
  fscanf(pFile,"%d",&Np);
  
  if(Nu != Np){
  	printf("Non-Matching number of elements - Probably wrong pair of files\n"); 
  	fclose(uFile); fclose(pFile);
  }

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  fgets(buffer,1024,uFile);
  fgets(buffer,1024,pFile);

  for(i=0;i<N;i+=1){
    //fgets(buffer,1024,uFile);
    printf("%s",buffer);
  	fscanf(uFile," (%f%f%f)");
    //sscanf(buffer,"(%f%f%f)",&(v[i].ux),&(v[i].uy),&(v[i].uz));
  }

  for(i=0;i<N;i+=1)
  	fscanf(pFile,"%f",&(v[i].p));

  printf("printing read values\n");
  for(i=0;i<N;i+=1)
  	printf("%d) %f %f %f %f\n",i,v[i].ux,v[i].uy,v[i].uz,v[i].p);

  free(uFile); free(pFile);

  return 0;
}
