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

#define DEBUG_MODE true
#define DEBUG_PRINT true
#define SUBTRACTION_MODE true

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

int readAxis(int Nx,int Ny,int Nz,int planeType,
             char *folder, double *X,double *Y);

int foamCalcScalar(int runType,int calcScalarMode,int Height,int Width,
                   int padWidth,double *X,double *Y,double *Xbuff,double *Ybuff,
                   double *uField,double *uSubtr,double *uBuff,
                   double *ux,double *uy,double *uxxx,double *uyyy,
                   double *uxxy,double *uxyy,double *gField,double *g2Field,
                   double v0y0,double *sField,double *sSubtr);

int main(int argc,char **argv){
  //double cutoff; int rCnet=0;
  int Width = 100, Height = 100, Depth,nFixVortex=5,NumCls=2048;
  int dataSize=12,runType=0,*label=NULL,**eqClass=NULL,nSkip=0;
  int *eqPop, *labelTag;
  int Nsnapshots=0,openFoamFile=0;
  int Nx,Ny,Nz,planeIndex,planeType,planeNum=0,pln[8128];
  int i,j,l,err,nCnect=0,n,nMax=1024,padWidth=2,calcScalarMode;//,k;
  double Gmin=1.,Gmax=20.,rmin=0.5,rmax=1.0,threshold=0.5;
  double xmin[2]={-9.,-9.},xmax[2]={9.,9.},x0[2],dx[2],xf[2],*background;
  double *parVortex=NULL,t,t0,dt;
  double *sField=NULL,*gField=NULL,*g2Field=NULL,*uField=NULL,*X,*Y,*wBkg;
  double *uBuff=NULL,*Xbuff=NULL,*Ybuff=NULL,*Xload=NULL,*Yload=NULL;
  double *Zload=NULL,*ux=NULL,*uy=NULL,*uxxy=NULL,*uxyy=NULL;
  double *uxxx=NULL,*uyyy=NULL,*vCatalog=NULL,*rCatalog=NULL,*uVort;
  double v0y0 = 0.00,*vortSndMomMatrix=NULL,*avgGradU=NULL,*sSubtr,*uSubtr;
  char folder[100+1],tag[100+1],filename[400+1],foamFolder[200+1],bkgFile[400+1];
  FILE *vortexFile,*dadosin;
  configVar cfg;
  openFoamIcoData *node=NULL; 
  ///////////////
  int iT=0,Tw=0,chunk;
  int iX=0,iY=0,iZ=0,Xw=0,Yw=0,Zw=0;
  char authtoken[400+1];
  char dataset[400+1],jhtdbFolder[400+1];
  ///////////////

  NumCls = 65536;

  dataSize = 6;//6;

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

  /// JHU data

  iX = cfg.jhtdb_Raw_iX; iY = cfg.jhtdb_Raw_iY; iZ = cfg.jhtdb_Raw_iZ;
  Xw = cfg.jhtdb_Raw_Xw; Yw = cfg.jhtdb_Raw_Yw; Zw = cfg.jhtdb_Raw_Zw;
  //strcpy(dataset,cfg.jhtdb_dataset);
  //strcpy(jhtdbFolder,cfg.jhtdb_folder); 
  //strcpy(authtoken,cfg.jhtdb_authToken);
  
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
    return 4;
    planeNum = 1;
    pln[0] = planeIndex;
  }

       if(planeType==0){ Height = Ny; Width = Nx; Depth = Nz; }
  else if(planeType==1){ Height = Ny; Width = Nz; Depth = Nx; }
  else if(planeType==2){ Height = Nz; Width = Nx; Depth = Ny; }
  else if(planeType==3){ Height = Yw; Width = Xw; Depth = 1 ; }
  else if(planeType==4){ Height = Yw; Width = Zw; Depth = 1 ; }
  else if(planeType==5){ Height = Zw; Width = Xw; Depth = 1 ; }
  else{ printf("error, non-recognized plane type\n"); return -15; }

  if(planeIndex<0)
    printf("Switching to plane number list\n");

  if(planeIndex>Depth)
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

  dbgPrint(10,0);

  /* Memory Allocation */

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
  fieldAlloc(wBkg,Height*Width,double);
  fieldAlloc(uBuff ,2*(Height+2*padWidth)*(Width+2*padWidth),double);
  fieldAlloc(X,Nx+1,double);
  fieldAlloc(Y,Ny+1,double);
  fieldAlloc(Xbuff,Width+2*padWidth,double);
  fieldAlloc(Ybuff,Height+2*padWidth,double);
  fieldAlloc(Xload,Nx+1,double);
  fieldAlloc(Yload,Ny+1,double);
  fieldAlloc(Zload,Nz+1,double);
  fieldAlloc(rCatalog,dataSize*nMax,double);
  fieldAlloc(vCatalog,4*nMax,double);
  fieldAlloc(uVort,2*nMax,double);
  fieldAlloc(vortSndMomMatrix,4*nMax,double);
  fieldAlloc(avgGradU,4*nMax,double);

  node = (openFoamIcoData*)malloc(Nx*Ny*Nz*sizeof(openFoamIcoData));
  if(node==NULL){
    printf("problems alocating openFoamIcoData node\n");
  }

  eqClass=(int**)malloc(NumCls*sizeof(int*));
  if(eqClass==NULL)
    return 1;
  for(i=0;i<NumCls;i+=1){
    eqClass[i]=(int*)malloc(NumCls*sizeof(int));
    if(eqClass[i]==NULL)
      return(i+2);
  }
  
  eqPop = (int*)malloc(NumCls*sizeof(int));
  if(eqPop==NULL)
    return 1;

  labelTag = (int*)malloc(NumCls*sizeof(int));
  if(labelTag==NULL)
    return 1;
  
  dbgPrint(13,0);

  if(DEBUG_MODE==true){
    printf("%d %d %d \n",Height,Width,nFixVortex);
    printf("%f %f %f %f %f %f\n",x0[0],x0[1],xf[0],xf[1],dx[0],dx[1]);
    printf("%f %f %f %f\n",xmin[0],xmin[1],xmax[0],xmax[1]);
    printf("%f %f %f %f\n",Gmin,Gmax,rmin,rmax);
  }

  dbgPrint(14,0);
  
  err=readAxis(Nx,Ny,Nz,planeType,folder,X,Y);
  if(DEBUG_MODE){
    for(i=0;i<Height;i+=1)
      printf("%lf\n",Y[i]);
    printf("\n");

    for(j=0;j<Width;j+=1)
      printf("%lf\n",X[j]);
    printf("\n");
  }

  dbgPrint(14,1);

  err = XtoXbuff(Width,X,Xbuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - X\n");

  err = XtoXbuff(Height,Y,Ybuff,padWidth);
  if(err!=0)
    printf("problem in XtoXbuff - Y\n");

  dbgPrint(14,2);

  sprintf(filename,"%s/vortices-%s.dat",folder,tag);
  vortexFile = fopen(filename,"w");

  dbgPrint(14,3);

  if(bkgFile[0]!='\0'){
    double Ux,Uy;
    FILE *dadosField;
    if(DEBUG_PRINT)
      printf("loading background file\n");
    dadosField=fopen(bkgFile,"r");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        fscanf(dadosField,"%lf %lf",&Ux,&Uy);
        background[2*(i*Width+j)+0] = Ux;
        background[2*(i*Width+j)+1] = Uy;
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
    }

    for(l=0;l<planeNum;l+=1){
      
      openFoamFile = 0;
      for(i=0;i<2*Height*Width;i+=1)
        uField[i]=0.;

      dbgPrint(15,0);

      if(DEBUG_PRINT)
        printf("plane =%d\n",pln[l]);

      if(planeType==0)      sprintf(filename,"%s/plane-z%d-%g.dat",folder,pln[l],t);
      else if(planeType==1) sprintf(filename,"%s/plane-x%d-%g.dat",folder,pln[l],t);
      else if(planeType==2) sprintf(filename,"%s/plane-y%d-%g.dat",folder,pln[l],t);
      else if(planeType==3) sprintf(filename,"%s/slice-(%d,%d,%d)-(%d,%d,%d)-%d.dat",folder,iX,iY,pln[l],Xw,Yw,Zw,n);
      else if(planeType==4) sprintf(filename,"%s/slice-(%d,%d,%d)-(%d,%d,%d)-%d.dat",folder,iX,pln[l],iZ,Xw,Yw,Zw,n);
      else if(planeType==5) sprintf(filename,"%s/slice-(%d,%d,%d)-(%d,%d,%d)-%d.dat",folder,pln[l],iY,iZ,Xw,Yw,Zw,n);
      
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
        
      if(openFoamFile!=0)
        continue;

      dbgPrint(15,3);

      for(i=0;i<Height;i+=1)
        for(j=0;j<Width;j+=1){
          uSubtr[2*(i*Width+j)+0]= uField[2*(i*Width+j)+0] 
                                 - background[2*(i*Width+j)+0];
          uSubtr[2*(i*Width+j)+1]= uField[2*(i*Width+j)+1]
                                 - background[2*(i*Width+j)+1];
        }
      
      
      err=foamCalcScalar(runType,calcScalarMode,Height,Width,padWidth,X,Y,
                         Xbuff,Ybuff,uField,uSubtr,uBuff,ux,uy,uxxx,uyyy,uxxy,
                         uxyy,gField,g2Field,v0y0,sField,sSubtr);
      if(err==-1)
        break;
      if(err!=0){
          printf("Error in calcScalarField - %d\n",err);
          return err;
        }
      
      dbgPrint(15,8);

      for(i=0;i<Height*Width;i+=1)
        label[i]=-1;

      err = floodFill(sField,Width,Height,NumCls,eqClass,eqPop,label);
      if(err!=0)
        printf("Problems in floodFill\n");

      err = renameLabels(Height,Width,NumCls,labelTag,label);
      if(err>0)
        nCnect=err;
      else
        printf("problems with renameLabels - %d\n",err);

      dbgPrint(15,9);
      
      if(n%100==0)
        printScalarFields(Height,Width,X,Y,n,folder,pln[l],sField,label);
      
      dbgPrint(15,10);

      err=extLambOseenParams(Height,Width,nCnect,X,Y,sField,gField,label,
                             vCatalog);
      if(err!=0){
        printf("problems in extLambOseenParams\n");
        return err;
      }

      dbgPrint(15,11);   
      
      err=extVortexVelocity(Height,Width,nCnect,X,Y,uField,sField,
                            gField,label,uVort);
      if(err!=0){
        printf("problems in extVortexVelocity\n");
        return err;
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
        rCatalog[dataSize*i+0] = vCatalog[4*i+0];
        rCatalog[dataSize*i+1] = vCatalog[4*i+1];
        rCatalog[dataSize*i+2] = vCatalog[4*i+2];
        rCatalog[dataSize*i+3] = vCatalog[4*i+3];
        // Added as to add vortex avg velocity
        rCatalog[dataSize*i+4] = uVort[2*i+0];
        rCatalog[dataSize*i+5] = uVort[2*i+1];
      }
  
      vortexAdaptiveQuickSort(rCatalog,nCnect,dataSize,&greaterAbsCirculation);
    
      dbgPrint(18,0);

      /* Preparing for printing */

      err=fprintVortex(vortexFile,n,dataSize,nCnect,rCatalog);
      if(err!=0){printf("problems in printing vortexfile\n"); return -6;}
    }
  } // End of Main loop

    
  printf("%d timesteps processed\n",n);

  fclose(vortexFile);
    
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
  if(avgGradU!=NULL) free(avgGradU);
  if(vortSndMomMatrix!=NULL) free(vortSndMomMatrix);
  if(uVort!=NULL) free(uVort);
  if(parVortex!=0) free(parVortex);
  if(vCatalog!=NULL) free(vCatalog);
  if(rCatalog!=NULL) free(rCatalog);

  dbgPrint(23,0);
 
  for(i=0;i<NumCls;i+=1)
    free(eqClass[i]);
  free(eqClass);

  dbgPrint(24,0);

  return 0;
}

int readAxis(int Nx,int Ny,int Nz,int planeType,
             char *folder, double *X,double *Y)
{
  int i;
  char filename[400+1];
  FILE *nFile;

  if(planeType==0 || planeType==3){
    sprintf(filename,"%s/Xaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Yaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }
  else if(planeType==1 || planeType==4){
    sprintf(filename,"%s/Zaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Yaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }
  else if(planeType==2 || planeType==5){
    sprintf(filename,"%s/Zaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Nx;i+=1)
      fscanf(nFile,"%lf",&(X[i]));
    fclose(nFile); nFile=NULL;

    sprintf(filename,"%s/Xaxis.dat",folder);
    nFile=fopen(filename,"r");
    for(i=0;i<Ny;i+=1)
      fscanf(nFile,"%lf",&(Y[i]));
    fclose(nFile); nFile=NULL;
  }

  return 0;
}

int foamCalcScalar(int runType,int calcScalarMode,int Height,int Width,
                   int padWidth,double *X,double *Y,double *Xbuff,double *Ybuff,
                   double *uField,double *uSubtr,double *uBuff,
                   double *ux,double *uy,double *uxxx,double *uyyy,
                   double *uxxy,double *uxyy,double *gField,double *g2Field,
                   double v0y0,double *sField,double *sSubtr)
{
  int err,i,j;
 
  if(calcScalarMode==0){
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sField);
    if(err!=0)
      return -1;
  }
  else if(calcScalarMode==2){
    dbgPrint(15,4);
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uSubtr,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sSubtr);
    if(err!=0)
      return -1;

    dbgPrint(15,5);
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sField);
    if(err!=0)
      return -1;
 
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        if( sField[i*Width+j] > sSubtr[i*Width+j])
          sField[i*Width+j] = sSubtr[i*Width+j];
        
        // WARNING : Why the hell are this velocity substitution by its subtraction is here?
        //uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
        //uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
      }
    dbgPrint(15,6);
  }
  else if(calcScalarMode==3){
    dbgPrint(15,4);
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sField);
    if(err!=0)
      return -1;

    dbgPrint(15,6);
    err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                        uSubtr,uBuff,ux,uy,uxxx,uyyy,uxxy,
                        uxyy,gField,g2Field,v0y0,sSubtr);
    if(err!=0)
      return -1;

    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        if( (sField[i*Width+j]==0) || (sSubtr[i*Width+j]==0))
          sField[i*Width+j] = 0.;

        // WARNING : Why the hell are this velocity substitution by its subtraction is here?
        //uField[2*(i*Width+j)+0]= uSubtr[2*(i*Width+j)+0];
        //uField[2*(i*Width+j)+1]= uSubtr[2*(i*Width+j)+1];
      }
  }
  else{
    printf("Not identified operation mode - %d\n",calcScalarMode);
    return -18;
  }

  return 0;
} 
