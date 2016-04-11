#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));\
                                  if(ptr==NULL){                         \
                                    printf("memory not allocked\n");     \
                                    return 1;                            \
                                  }                                      \
                                  else{                                  \
                                    for(i=0;i<(size);i+=1)               \
                                      ptr[i]=(type) 0;                   \
                                  }                                      \

int main(int argc,char **argv){
  const int Npre=20;
  int i,n,N,Nsteps,Nskip,Nu,Np;
  double t0,dt,t;
  double *uAvg,*pAvg,u,v,w,p;
  char buffer[1024],filename[1024];
  FILE *uFile,*pFile,*uSubFile,*pSubFile;
  
  Nsteps = 10;

  t0 = 10.0075;//atod(argv[2]);
  dt =  0.0120;//atod(argv[3]);

  N=256*192*192;
  printf("N=%d\n",N);

  fieldAlloc(uAvg,3*N,double);
  fieldAlloc(pAvg,N,double);
  
  Nskip = 0;

  for(n=0;n<Nsteps;n+=1){
    t=t0+((double)n)*dt;

    printf("t=%lf\n",t);
    
    sprintf(filename,"%s/%g/U",argv[1],t);
    uFile = fopen(filename,"r");
    
    sprintf(filename,"%s/%g/p",argv[1],t);
    pFile = fopen(filename,"r");

    if(uFile==NULL || pFile==NULL){ 
      Nskip +=1;
      printf("problems opening uFile or pFile- %d\n",n);
      continue;
    }

    for(i=0;i<Npre;i++)
      fgets(buffer,1024,uFile);
    fscanf(uFile,"%d",&Nu);

    for(i=0;i<Npre;i++)
      fgets(buffer,1024,pFile);
    fscanf(pFile,"%d",&Np);
  
    if(Nu != Np){
      fclose(uFile); fclose(pFile);
      return -10;
    }

    fgets(buffer,1024,uFile);
    fgets(buffer,1024,pFile);

    fgets(buffer,1024,uFile);
    fgets(buffer,1024,pFile);
    
    for(i=0;i<N;i+=1){
      fscanf(pFile,"%lf",&p);
      fscanf(uFile," (%lf%lf%lf)",&u,&v,&w);

      pAvg[i] += p;
      uAvg[3*i+0] += u;
      uAvg[3*i+1] += v;
      uAvg[3*i+2] += w;
    }
    
    if(pFile!=NULL){fclose(pFile);pFile=NULL;}
    if(uFile!=NULL){fclose(uFile);uFile=NULL;}
  }

  for(i=0;i<N;i+=1){
  	pAvg[i] /= Nsteps-Nskip;
    uAvg[3*i+0] /= Nsteps-Nskip;
    uAvg[3*i+1] /= Nsteps-Nskip;
    uAvg[3*i+2] /= Nsteps-Nskip;
  }

  sprintf(filename,"avgU.dat");
  uFile = fopen(filename,"w");
    
  sprintf(filename,"avgp.dat");
  pFile = fopen(filename,"w");

  fprintf(uFile,
  	      "/*--------------------------------*- C++ -*----------------------------------*\\\n"
          "| =========                 |                                                 | \n"
          "| \\\\      /  F ield         | OpenFOAM Extend Project: Open source CFD        | \n"
          "|  \\\\    /   O peration     | Version:  1.6-ext                               | \n"
          "|   \\\\  /    A nd           | Web:      www.extend-project.de                 | \n"
          "|    \\\\/     M anipulation  |                                                 | \n"
          "\\*---------------------------------------------------------------------------*/\n"
          "FoamFile\n"
          "{\n"
          "    version     2.0;\n"
          "    format      ascii;\n"
          "    class       volVectorField;\n"
          "    location    \"0.0\";\n"
          "    object      U;\n"
          "}\n"
          "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
          "\n"
          "dimensions      [0 1 -1 0 0 0 0];\n"
          "internalField   nonuniform List<vector>\n"
          "%d\n"
          "(\n",N);

  for(i=0;i<N;i+=1)
  	fprintf(uFile,"(%lf %lf %lf)\n",uAvg[3*i+0],uAvg[3*i+1],uAvg[3*i+2]);
  
  fprintf(uFile,
  	")\n"
    ";\n"
    "    }\n"
    "}\n"
    "\n"
    "\n"
    "// ************************************************************************* //\n"
    "\n"
    "\n"
  );

  fprintf(pFile,
  	      "/*--------------------------------*- C++ -*----------------------------------*\\\n"
          "| =========                 |                                                 | \n"
          "| \\\\      /  F ield         | OpenFOAM Extend Project: Open source CFD        | \n"
          "|  \\\\    /   O peration     | Version:  1.6-ext                               | \n"
          "|   \\\\  /    A nd           | Web:      www.extend-project.de                 | \n"
          "|    \\\\/     M anipulation  |                                                 | \n"
          "\\*---------------------------------------------------------------------------*/\n"
          "FoamFile\n"
          "{\n"
          "    version     2.0;\n"
          "    format      ascii;\n"
          "    class       volScalarField;\n"
          "    location    \"0.0\";\n"
          "    object      p;\n"
          "}\n"
          "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
          "\n"
          "dimensions      [0 2 -2 0 0 0 0];\n"
          "internalField   nonuniform List<scalar>\n"
          "%d\n"
          "(\n",N);

  for(i=0;i<N;i+=1)
  	fprintf(pFile,"%lf\n",pAvg[i]);
  
  fprintf(pFile,
          ")"
          ";");

  fclose(pFile); fclose(uFile);

  printf("adding subtracted field\n");
  for(n=0;n<Nsteps;n+=1){
    t=t0+((double)n)*dt;

    printf("t=%lf\n",t);
    
    sprintf(filename,"%s/%g/U",argv[1],t);
    uFile = fopen(filename,"r");
    
    sprintf(filename,"%s/%g/p",argv[1],t);
    pFile = fopen(filename,"r");

    sprintf(filename,"%s/%g/Usub",argv[1],t);
    uSubFile = fopen(filename,"w");
    
    sprintf(filename,"%s/%g/psub",argv[1],t);
    pSubFile = fopen(filename,"w");

    if(uFile==NULL || pFile==NULL || uSubFile==NULL || pSubFile==NULL){ 
      Nskip +=1;
      printf("problems opening uFile or pFile or uSubFile or pSubFile- %d\n",n);
      continue;
    }

    for(i=0;i<Npre;i++)
      fgets(buffer,1024,uFile);
    fscanf(uFile,"%d",&Nu);

    for(i=0;i<Npre;i++)
      fgets(buffer,1024,pFile);
    fscanf(pFile,"%d",&Np);
  
    if(Nu != Np){
      fclose(uFile); fclose(pFile);
      return -10;
    }

    fgets(buffer,1024,uFile);
    fgets(buffer,1024,pFile);

    fgets(buffer,1024,uFile);
    fgets(buffer,1024,pFile);
    
    fprintf(uSubFile,
  	        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
            "| =========                 |                                                 | \n"
            "| \\\\      /  F ield         | OpenFOAM Extend Project: Open source CFD        | \n"
            "|  \\\\    /   O peration     | Version:  1.6-ext                               | \n"
            "|   \\\\  /    A nd           | Web:      www.extend-project.de                 | \n"
            "|    \\\\/     M anipulation  |                                                 | \n"
            "\\*---------------------------------------------------------------------------*/\n"
            "FoamFile\n"
            "{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       volVectorField;\n"
            "    location    \"%f\";\n"
            "    object      U;\n"
            "}\n"
            "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
            "\n"
            "dimensions      [0 1 -1 0 0 0 0];\n"
            "internalField   nonuniform List<vector>\n"
            "%d\n"
            "(\n",t,N);
    
    fprintf(pSubFile,
  	        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
            "| =========                 |                                                 | \n"
            "| \\\\      /  F ield         | OpenFOAM Extend Project: Open source CFD        | \n"
            "|  \\\\    /   O peration     | Version:  1.6-ext                               | \n"
            "|   \\\\  /    A nd           | Web:      www.extend-project.de                 | \n"
            "|    \\\\/     M anipulation  |                                                 | \n"
            "\\*---------------------------------------------------------------------------*/\n"
            "FoamFile\n"
            "{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       volScalarField;\n"
            "    location    \"%f\";\n"
            "    object      p;\n"
            "}\n"
            "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
            "\n"
            "dimensions      [0 2 -2 0 0 0 0];\n"
            "internalField   nonuniform List<scalar>\n"
            "%d\n"
            "(\n",t,N);

    for(i=0;i<N;i+=1){
      fscanf(pFile,"%lf",&p);
      fscanf(uFile," (%lf%lf%lf)",&u,&v,&w);
      
      fprintf(uSubFile,"(%lf %lf %lf)\n",u-uAvg[3*i+0],v-uAvg[3*i+1],w-uAvg[3*i+2]);
      fprintf(pSubFile,"%lf\n",p-pAvg[i]);
    }


    fprintf(uSubFile,
            ")\n"
            ";\n"
            "    }\n"
            "}\n"
            "\n"
            "\n"
            "// ************************************************************************* //\n"
            "\n" 
            "\n"
    );

    fprintf(pSubFile,
            ")"
            ";");
    
    if(pFile!=NULL){fclose(pFile);pFile=NULL;}
    if(uFile!=NULL){fclose(uFile);uFile=NULL;}
    if(pSubFile!=NULL){fclose(pSubFile);pFile=NULL;}
    if(uSubFile!=NULL){fclose(uSubFile);uFile=NULL;}
  }

  return 0;
}