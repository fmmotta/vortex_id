#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_histogram.h>
#include "ini.h"
#include "mt64.h"
#include "vortexGen.h"
#include "floodFill.h"
#include "stencilExtended.h"
#include "preprocessing.h"
#include "lambdaInit.h"
#include "vortexExtraction.h"
#include "vortexExtractionExtend.h"
#include "inputManager.h"
#include "essayHandler.h"

#define DEBUG_MODE false
#define DEBUG_PRINT false
#define SUBTRACTION_MODE false

#define dbgPrint(num,num2) if(DEBUG_PRINT) printf("check point - %d-%d\n",(num),(num2))

#define fieldAlloc(ptr,size,type) ptr=(type*)malloc((size)*sizeof(type));\
                                  if(ptr==NULL){                         \
                                    printf("memory not allocked\n");     \
                                    return 1;                            \
                                  }                                      \
                                  else{                                  \
                                    for(i=0;i<(size);i+=1)               \
                                      ptr[i]=(type) 0;                   \
                                  }                                      \

#define sfopen(file,name,mode) file = fopen(name,mode);                     \
                               if(file==NULL)                               \
                                printf("problems opening file: %s\n",name); 

int main(int argc,char **argv){
  //double cutoff; int rCnet=0;
  int Width = 100, Height = 100, Depth, nVortex=5,nFixVortex=5;
  int dataSize=12,runType=0,*label=NULL,**eqClass=NULL,nSkip=0;
  int hNG=50,hNRc=53,hNa=40,hNb=40,hNN=10,Nsnapshots=0,openFoamFile=0;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum=0,pln[8128];
  int i,j,l,err,nCnect=0,n,nMax=1024,padWidth=2,calcScalarMode;//,k;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2],*background;
  double *parVortex=NULL,t,t0,dt,*avgU=NULL,*avgU2=NULL,*avgW=NULL,*avgW2=NULL;
  double *sField=NULL,*gField=NULL,*g2Field=NULL,*uField=NULL,*X,*Y,*wBkg;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*Xload=NULL,*Yload=NULL;
  double *Zload=NULL,*ux=NULL,*uy=NULL,*uxxy=NULL,*uxyy=NULL,*sSubtr=NULL;
  double *uxxx=NULL,*uyyy=NULL,*vCatalog=NULL,*rCatalog=NULL,*uSubtr=NULL;
  double v0y0 = 0.00,*vortSndMomMatrix=NULL,*avgGradU=NULL;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1],bkgFile[400+1];
  FILE *dadosout,*uFile,*pFile,*nFile,*vortexFile,*totalVortices,*dadosin;
  configVar cfg;
  openFoamIcoData *node=NULL; 

  dataSize = 12;

  if(argc!=2){
    printf("Incorrect Number of Arguments - Need exactly "
           "the configuration file\n");
    return -1;
  }

  err=initConfig(&cfg);
  
  if (ini_parse(argv[1], vortexIdHandler, &cfg) < 0){
    printf("Can't load .ini file\n");
    return 1;
  }
  
  dbgPrint(2,0);

  if(DEBUG_MODE==true){
    err=printConfig(&cfg);
    if(err!=0)
      return err;
  }
  
  if(cfg.dim!=2){
    printf("Dimension is not 2 - Can't Follow\n");
    return 2;
  }
  
  /**** OpenFOAM parameters setting *****/
  
  Nx = cfg.Nx;
  Ny = cfg.Ny;
  Nz = cfg.Nz;
  t0 = cfg.t0;
  dt = cfg.dt;
  planeIndex = cfg.pIndex;
  planeType  = cfg.pType;
  Nsnapshots = cfg.Nsnapshots;
  strcpy(foamFolder,cfg.FOAMfolder);
  
  if(planeIndex<0){
    if(cfg.planeNum>0){
      planeNum=cfg.planeNum;
      for(i=0;i<planeNum;i+=1)
        pln[i]=cfg.pln[i];
    }
    else{
      printf("Wrongly written configuration file, specify number of slices\n");
      return 1;
    }
  }else{
  	printf("This should not happen\n");
    planeNum = 1;
    pln[0] = planeIndex;
  }

  if(planeType==0){
    Height = Ny;
    Width  = Nx;
    Depth  = Nz;
  }
  else if(planeType==1){
    Height = Ny;
    Width  = Nz;
    Depth  = Nx; 
  }
  else if(planeType==2){
    Height = Nz;
    Width  = Nx;
    Depth  = Ny; 
  }
  else{
    printf("error, non-recognized plane type\n");
    return -15;
  }

  if(planeIndex<0)
    printf("Switching to plane number list\n");

  if(planeIndex>=Depth)
    printf("Out of bounds plane\n");

  if(cfg.Nx == 0 || cfg.Ny == 0 || cfg.Nz == 0){
    printf("error, incompatible dimension sizes\n");
    return -16;
  }

  /* Loading Configuration -- I need something more concise */

  dbgPrint(2,0);

  strcpy(folder,cfg.folder);
  strcpy(tag,cfg.tag);

  if(cfg.bkgFile!=NULL)
    strcpy(bkgFile,cfg.bkgFile);
  else
    bkgFile[0]='\0';

  err=mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(err!=0 && err!=-1){
    printf("error creating directory - %d\n",err);
    return err;
  }
  
  runType = cfg.runType;
  calcScalarMode = cfg.calcMode;
  
  dbgPrint(4,1);

  /**********************************/

  dbgPrint(7,0);
  
  v0y0    = cfg.v0y0;
  //cutoff  = cfg.cutoff;

  if(runType==0)
    threshold = cfg.swThresh;
  else if (runType==2)
    threshold = cfg.sndSwThresh;
  else
    threshold = -1.;

  if(DEBUG_MODE)
    printf("threshold is %f\n",threshold);
  
  dbgPrint(8,0);

  err=freeConfig(&cfg);if(err!=0) return err;

  /* End Loading Configuration */
  /* Memory Allocation */

  dbgPrint(10,0);

  fieldAlloc(label,Height*Width,int);
  fieldAlloc(sField ,Height*Width,double);
  fieldAlloc(sSubtr ,Height*Width,double);
  fieldAlloc(gField ,4*Height*Width,double);
  fieldAlloc(g2Field,4*Height*Width,double);
  fieldAlloc(uField,2*Height*Width,double);
  fieldAlloc(uSubtr,2*Height*Width,double);
  fieldAlloc(background,2*Height*Width,double);
  fieldAlloc(  ux  ,2*Height*Width,double);
  fieldAlloc(  uy  ,2*Height*Width,double);
  fieldAlloc( uxxy ,2*Height*Width,double);
  fieldAlloc( uxyy ,2*Height*Width,double);
  fieldAlloc( uxxx ,2*Height*Width,double);
  fieldAlloc( uyyy ,2*Height*Width,double);
  fieldAlloc(avgU,2*Height*Width,double);
  fieldAlloc(avgU2,2*Height*Width,double);
  fieldAlloc(avgW,Height*Width,double);
  fieldAlloc(avgW2,Height*Width,double);
  fieldAlloc(wBkg,Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);
  fieldAlloc(X,Nx+1,double);
  fieldAlloc(Y,Ny+1,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);
  fieldAlloc(Xload,Nx+1,double);
  fieldAlloc(Yload,Ny+1,double);
  fieldAlloc(Zload,Nz+1,double);
  
  dbgPrint(5,0);

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }

  vCatalog = (double*)malloc(4*nMax*sizeof(double));
  if(vCatalog==NULL){
    printf("memory not allocked\n");
    return 3;
  }
  for(i=0;i<4*nMax;i+=1)
    vCatalog[i]=-1.;
  
  rCatalog = (double*)malloc(dataSize*nMax*sizeof(double));
  if(rCatalog==NULL){
    printf("memory not allocked\n");
    return 4;
  }
  for(i=0;i<4*nMax;i+=1)
    rCatalog[i]=-1.;
  
  vortSndMomMatrix = (double*)malloc(4*nMax*sizeof(double));
  if(vortSndMomMatrix==NULL){
    printf("memory not allocked\n");
    return 3;
  }
  for(i=0;i<4*nMax;i+=1)
    vortSndMomMatrix[i]=-1.;

  avgGradU = (double*)malloc(4*nMax*sizeof(double));
  if(avgGradU==NULL){
    printf("memory not allocked\n");
    return 3;
  }
  for(i=0;i<4*nMax;i+=1)
    avgGradU[i]=-0.;

  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("not enough memory for openFoamIcoData\n");
    return 1;
  }

  dbgPrint(13,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  dbgPrint(14,0);

  if(planeIndex>=0){
    sprintf(filename,"%s/constant/polyMesh/points",foamFolder);
    sfopen(nFile,filename,"r");
    err=loadAxis(nFile,Nx,Ny,Nz,Xload,Yload,Zload);
    if(err!=0)
      return err;
    fclose(nFile);

    if(planeType==0){
      for(i=0;i<Height;i+=1)
        Y[i] = (Yload[i]+Yload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Xload[j]+Xload[j+1])/2.;
    }
    else if(planeType==1){
      for(i=0;i<Height;i+=1)
        Y[i] = (Yload[i]+Yload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Zload[j]+Zload[j+1])/2.; 
    }
    else if(planeType==2){
      for(i=0;i<Height;i+=1)
        Y[i] = (Zload[i]+Zload[i+1])/2.;
      for(j=0;j<Width;j+=1)
        X[j] = (Xload[j]+Xload[j+1])/2.; 
    }
    else 
      printf("non-identified plane type\n");
  }
  else{
    if(planeType==0){
      sprintf(filename,"%s/Xaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Nx;i+=1)
        fscanf(nFile,"%lf",&(X[i]));
      fclose(nFile); nFile=NULL;

      sprintf(filename,"%s/Yaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Ny;i+=1)
        fscanf(nFile,"%lf",&(Y[i]));
      fclose(nFile); nFile=NULL;
    }
    else if(planeType==1){
      sprintf(filename,"%s/Yaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Nx;i+=1)
        fscanf(nFile,"%lf",&(X[i]));
      fclose(nFile); nFile=NULL;

      sprintf(filename,"%s/Zaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Ny;i+=1)
        fscanf(nFile,"%lf",&(Y[i]));
      fclose(nFile); nFile=NULL;
    }
    else if(planeType==2){
      sprintf(filename,"%s/Zaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Nx;i+=1)
        fscanf(nFile,"%lf",&(X[i]));
      fclose(nFile); nFile=NULL;

      sprintf(filename,"%s/Xaxis.dat",folder);
      sfopen(nFile,filename,"r");
      for(i=0;i<Ny;i+=1)
        fscanf(nFile,"%lf",&(Y[i]));
      fclose(nFile); nFile=NULL;
    }
  }

  dbgPrint(14,1);

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  dbgPrint(14,2);

  sprintf(filename,"%s/vortices.dat",folder);
  sfopen(vortexFile,filename,"w");
  
  sprintf(filename,"%s/totalVortices.dat",folder);
  sfopen(totalVortices,filename,"w");

  dbgPrint(14,3);

  for(i=0;i<2*Height*Width;i+=1)
    avgU[i]=0.;
  for(i=0;i<2*Height*Width;i+=1)
    avgU2[i]=0.;
  for(i=0;i<Height*Width;i+=1)
    avgW[i]=0.;
  for(i=0;i<Height*Width;i+=1)
    avgW2[i]=0.;
  for(i=0;i<2*Height*Width;i+=1)
    background[i]=0.;

  dbgPrint(14,4);

  if(bkgFile[0]!='\0'){
    double x,y,Ux,Uy,sigmaUx,sigmaUy;
    double omega,strain,gamma,beta;
    FILE *dadosField;
    if(DEBUG_PRINT)
      printf("loading background file\n");
    sfopen(dadosField,bkgFile,"r");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        fscanf(dadosField,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&Ux,&Uy,
                                                          &sigmaUx,&sigmaUy,
                                                          &omega,&gamma,
                                                          &beta,&strain);
        background[2*(i*Width+j)+0] = Ux;
        background[2*(i*Width+j)+1] = Uy;
        wBkg[i*Width+j] = omega;
      }

    fclose(dadosField);
  }

  dbgPrint(14,5);

  nSkip=0;

  for(n=0;n<Nsnapshots;n+=1){

    t=t0+((double)n)*dt;
    printf("%d timesteps processed\n",n);
    if(n%10 == 0){      
      fflush(vortexFile);
      fflush(totalVortices);
    }

    for(l=0;l<planeNum;l+=1){
      
      openFoamFile = 0;
      for(i=0;i<2*Height*Width;i+=1)
        uField[i]=0.;

      dbgPrint(15,0);

      if(planeIndex>0){

        //printf("planeIndex=%d\n",planeIndex);
        
        sprintf(filename,"%s/%g/U",foamFolder,t);
        sfopen(uFile,filename,"r");
        if(uFile==NULL) printf("problems opening uFile - %d\n",n);

        sprintf(filename,"%s/%g/p",foamFolder,t);
        sfopen(pFile,filename,"r");
        if(pFile==NULL) printf("problems opening pFile - %d\n",n);

        if(uFile == NULL || pFile == NULL){
          openFoamFile = 1;
          nSkip+=1;
          printf("Failed time = %g\n",t);
          break;
        }

        dbgPrint(15,1);
      
        err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
        if(err!=0) printf("Problems with loadFields\n");
      
        fclose(pFile); fclose(uFile);
         
        dbgPrint(15,2);
      
        err=sliceFoamField(Height,Width,planeType,pln[l],Nx,Ny,Nz,node,uField);
        if(err!=0)
          printf("Problem slicing openFoamIcoData\n");
      }
      else{
        if(DEBUG_PRINT)
          printf("Operating with multiple slices\n");
        
        if(DEBUG_PRINT)
          printf("plane =%d\n",pln[l]);

        if(planeType==0)      sprintf(filename,"%s/plane-z%d-%.4f.dat",folder,pln[l],t);
        else if(planeType==1) sprintf(filename,"%s/plane-x%d-%.4f.dat",folder,pln[l],t);
        else if(planeType==2) sprintf(filename,"%s/plane-y%d-%.4f.dat",folder,pln[l],t);

        dadosin=fopen(filename,"r");
        if(dadosin==NULL){
          openFoamFile=1;
          nSkip+=1;
          printf("Failed to open slice - %s\n",filename);
          break;
        }

        dbgPrint(15,31);
        
        for(i=0;i<Height;i+=1)
          for(j=0;j<Width;j+=1)
            fscanf(dadosin,"%lf%lf",&(uField[2*(i*Width+j)+0]),&(uField[2*(i*Width+j)+1]));
        fclose(dadosin);
        
        dbgPrint(15,32);
      }

      if(openFoamFile!=0)
        continue;

      dbgPrint(15,3);
      
      if(calcScalarMode==0){
        err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                            uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                            uxyy,gField,g2Field,v0y0,sField);
        if(err!=0)
          break;
      }
      else if(calcScalarMode==2){

        for(i=0;i<Height;i+=1)
          for(j=0;j<Width;j+=1){
            uSubtr[2*(i*Width+j)+0]= uField[2*(i*Width+j)+0] 
                                   - background[2*(i*Width+j)+0];
            uSubtr[2*(i*Width+j)+1]= uField[2*(i*Width+j)+1]
                                   - background[2*(i*Width+j)+1];
          }
           
        dbgPrint(15,4);
        err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                            uSubtr,uBuff,ux,uy,uxxx,uyyy,uxxy,
                            uxyy,gField,g2Field,v0y0,sSubtr);
        if(err!=0)
          break;

        //for(i=0;i<Height;i+=1)
        //  for(j=0;j<Width;j+=1){
        //    uField[2*(i*Width+j)+0]= 0.;
        //    uField[2*(i*Width+j)+1]= 0.;
        //  }

        dbgPrint(15,5);
        err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                            uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                            uxyy,gField,g2Field,v0y0,sField);
        if(err!=0)
          break;

        for(i=0;i<Height;i+=1)
          for(j=0;j<Width;j+=1){
            if( sField[i*Width+j] > sSubtr[i*Width+j])
              sField[i*Width+j] = sSubtr[i*Width+j];

            uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
            uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
          }
        dbgPrint(15,6);
      }
      else if(calcScalarMode==3){

        //for(i=0;i<Height;i+=1)
        //  for(j=0;j<Width;j+=1){
        //    uField[2*(i*Width+j)+0]= 0.;
        //    uField[2*(i*Width+j)+1]= 0.;
        //  }

        dbgPrint(15,4);
        err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                            uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                            uxyy,gField,g2Field,v0y0,sField);
        if(err!=0)
          break;
        
        dbgPrint(15,5);
        for(i=0;i<Height;i+=1)
          for(j=0;j<Width;j+=1){
            uSubtr[2*(i*Width+j)+0]= uField[2*(i*Width+j)+0] 
                                   - background[2*(i*Width+j)+0];
            uSubtr[2*(i*Width+j)+1]= uField[2*(i*Width+j)+1]
                                   - background[2*(i*Width+j)+1];
          }

        dbgPrint(15,6);
        err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                            uSubtr,uBuff,ux,uy,uxxx,uyyy,uxxy,
                            uxyy,gField,g2Field,v0y0,sSubtr);
        if(err!=0)
          break;

        for(i=0;i<Height;i+=1)
          for(j=0;j<Width;j+=1){
            //if( sField[i*Width+j] > sSubtr[i*Width+j])
            //  sField[i*Width+j] = sSubtr[i*Width+j];
            if( (sField[i*Width+j]==0) || (sSubtr[i*Width+j]==0))
              sField[i*Width+j] = 0.;

            uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
            uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
          }
      }
      else{
        printf("Not identified operation mode - %d\n",calcScalarMode);
        return -18;
      }
      
      if(err!=0){
          printf("Error in calcScalarField - %d\n",err);
          return err;
        }
      
      dbgPrint(15,7);

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          avgU[2*(i*Width+j)+0] += uField[2*(i*Width+j)+0];
          avgU[2*(i*Width+j)+1] += uField[2*(i*Width+j)+1];

          avgU2[2*(i*Width+j)+0] += uField[2*(i*Width+j)+0]*uField[2*(i*Width+j)+0];
          avgU2[2*(i*Width+j)+1] += uField[2*(i*Width+j)+1]*uField[2*(i*Width+j)+1];

          avgW[i*Width+j]  += gField[4*(i*Width+j)+2]-gField[4*(i*Width+j)+1];
          
          avgW2[i*Width+j] += (gField[4*(i*Width+j)+2]-gField[4*(i*Width+j)+1])*
                              (gField[4*(i*Width+j)+2]-gField[4*(i*Width+j)+1]);
        }

      dbgPrint(15,8);

      for(i=0;i<Height*Width;i+=1)
        label[i]=-1;

      err = floodFill(sField,Width,Height,eqClass,label);
      if(err!=0)
        printf("Problems in floodFill\n");

      err = renameLabels(Height,Width,label);
      if(err>0)
        nCnect=err;
      else
        printf("problems with renameLabels - %d\n",err);

      dbgPrint(15,9);
      
      if(n%100==0)
        printScalarFields(Height,Width,X,Y,n,folder,pln[l],sField,label);
      
      dbgPrint(15,10);

      //err=extract012Momentsw2(Height,Width,nCnect,X,Y,sField,gField,
      //                        label,vCatalog,vortSndMomMatrix,avgGradU);

      err=extract012Momentsw2(Height,Width,nCnect,X,Y,sField,gField,label,
                              vCatalog,vortSndMomMatrix,avgGradU);
      if(err!=0){
        printf("problems in extract012Momentsw2\n");
        return err;
      }

      if(SUBTRACTION_MODE){
        double bkgG[nCnect];
        
        for(i=0;i<nCnect;i+=1)
          bkgG[i] = 0.;
        
        err = extractAvgBkgVort(Height,Width,X,Y,nCnect,label,wBkg,bkgG);
        for(i=0;i<nCnect;i+=1)
          vCatalog[4*i+0] = vCatalog[4*i+0]-bkgG[i];
      }

      dbgPrint(16,0);

      for(i=0;i<nCnect;i+=1){
        if(runType==0){
          vCatalog[4*i+0]=1.397948086*vCatalog[4*i+0];
          vCatalog[4*i+1]= (1./1.12091)*vCatalog[4*i+1];
        }
        else if(runType==1){
          vCatalog[4*i+0]= 2.541494083*vCatalog[4*i+0];
          vCatalog[4*i+1]= (sqrt(2.))*vCatalog[4*i+1]; 
        }
      }

      for(i=0;i<nCnect;i+=1){
        rCatalog[dataSize*i+0]  = vCatalog[4*i+0];
        rCatalog[dataSize*i+1]  = vCatalog[4*i+1];
        rCatalog[dataSize*i+2]  = vCatalog[4*i+2];
        rCatalog[dataSize*i+3]  = vCatalog[4*i+3];
        rCatalog[dataSize*i+4]  = vortSndMomMatrix[4*i+0];
        rCatalog[dataSize*i+5]  = vortSndMomMatrix[4*i+1];
        rCatalog[dataSize*i+6]  = vortSndMomMatrix[4*i+2];
        rCatalog[dataSize*i+7]  = vortSndMomMatrix[4*i+3];
        rCatalog[dataSize*i+8]  = avgGradU[4*i+0];
        rCatalog[dataSize*i+9]  = avgGradU[4*i+1];
        rCatalog[dataSize*i+10] = avgGradU[4*i+2];
        rCatalog[dataSize*i+11] = avgGradU[4*i+3];
      }
  
      vortexAdaptiveQuickSort(rCatalog,nCnect,dataSize,&greaterAbsCirculation);
    
      dbgPrint(18,0);

      /* Preparing for printing */

      printVorticesAndMoments(Height,Width,X,Y,t,n,folder,pln[l],nCnect,vCatalog,rCatalog);

      err=fprintSafeVortexMoments(totalVortices,n,dataSize,nCnect,rCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems vortexSafeMoments Total\n"); return -6;}
    
      err=fprintSafeVortex(vortexFile,n,nCnect,vCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems in printing vortexfile\n"); return -6;}
    }
  } // End of Main loop

    
  printf("%d timesteps processed\n",n);

  fclose(vortexFile);
  fclose(totalVortices);
    
  dbgPrint(20,0);

  sprintf(filename,"gnuplot_script.gnu");
  err=writeGnuplotScript(filename,folder,tag,0,nVortex);
  if(err!=0){printf("Error printing gnuplot script\n");return err;}

  dbgPrint(22,0);

  if(X!=NULL) free(X);
  if(Y!=NULL) free(Y);
  if(Xbuff!=NULL) free(Xbuff);
  if(Ybuff!=NULL) free(Ybuff);
  if(Xload!=NULL) free(Xload);
  if(Yload!=NULL) free(Yload);
  if(Zload!=NULL) free(Zload);
  if(node !=NULL) free(node);
  if(label!=NULL) free(label);
  if(uField!=NULL) free(uField);
  if(uSubtr!=NULL) free(uSubtr);
  if(background!=NULL) free(background);
  if(sField!=NULL) free(sField); 
  if(sSubtr!=NULL) free(sSubtr);
  if(gField!=NULL)  free(gField);
  if(g2Field!=NULL) free(g2Field);
  if(uBuff!=NULL) free(uBuff);
  if(ux!=NULL) free(ux);
  if(uy!=NULL) free(uy);
  if(uxxx!=NULL) free(uxxx);
  if(uyyy!=NULL) free(uyyy);
  if(uxxy!=NULL) free(uxxy);
  if(uxyy!=NULL) free(uxyy);
  if(wBkg!=NULL) free(wBkg);
  if(avgU!=NULL) free(avgU);
  if(avgW!=NULL) free(avgW);
  if(avgU2!=NULL) free(avgU2);
  if(avgW2!=NULL) free(avgW2);
  if(avgGradU!=NULL) free(avgGradU);
  if(vortSndMomMatrix!=NULL) free(vortSndMomMatrix);
  if(parVortex!=0) free(parVortex);
  if(vCatalog!=NULL) free(vCatalog);
  if(rCatalog!=NULL) free(rCatalog);

  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  dbgPrint(24,0);

  return 0;
}