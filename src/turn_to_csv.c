#include <stdio.h>

int main(int argc,char **argv){
  int i,N;
  double G,r,x,y,Ux,Uy;
  FILE *infile,*outfile;

  if(argc!=4)
  	return 1;

  infile  = fopen(argv[0],"r");
  outfile = fopen(argv[1],"w");
  N = atoi(argv[2])
  for(i=0;i<N;i+=1){
  	fscanf(infile,"%lf%lf%lf%lf%lf%lf",&G,&r,&x,&y,&Ux,&Uy);
  	fprintf(outfile,"%lf,%lf,%lf,%lf,%lf,%lf\n",G,r,x,y,Ux,Uy);
  }

  return 0;
}