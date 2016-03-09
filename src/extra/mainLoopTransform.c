#define safeFopen(handle,mode,format,complement) sprintf(filename,format,complement);      \
                                                 handle = fopen(filename,mode);            \
                                                 if(handle == NULL)                        \
                                                   printf("Could not open: %s\n",filename);\

//---------------------------------------------------------------------

int sliceFoamField(int Height,int Width,int planeType,int index,
                   int Nx,int Ny,int Nz,openFoamIcoData *node,double *uField)
{
  int i,j,k;

  if(Height<=0 || Width<=0 || index <=0 || node == NULL || uField == NULL)
    return 1;

  if(planeType==0){
    if(DEBUG_PRINT)
      printf("XY plane\n");
    k=index;
    for(j=0;j<Height;j+=1)
      for(i=0;i<Width;i+=1){
        uField[2*(j*Width+i)+0] = node[id(i,j,k)].u;
        uField[2*(j*Width+i)+1] = node[id(i,j,k)].v;
      }
  }
  else if(planeType==1){
    if(DEBUG_PRINT)
      printf("ZY plane\n");
    i=index;
    for(j=0;j<Height;j+=1)
      for(k=0;k<Width;k+=1){
        uField[2*(j*Width+k)+0] = node[id(i,j,k)].w;
        uField[2*(j*Width+k)+1] = node[id(i,j,k)].v;
      }
  }
  else if(planeType==2){
    if(DEBUG_PRINT)
      printf("ZX plane\n");
    j=index;
    for(k=0;k<Height;k+=1)
      for(i=0;i<Width;i+=1){
        uField[2*(k*Width+i)+0] = node[id(i,j,k)].w;
        uField[2*(k*Width+i)+1] = node[id(i,j,k)].u;
      }
  }
  return 0;
}

//---------------------------------------------------------------------

void printScalarFields(int Height,int Width, double *X,double *Y,int n,
	                   char *folder,int index,double *sField,int *label)
{
  char filename[512+1];
  FILE *dadosout;

  sprintf(filename,"%s/sField-p%3d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUsfield(dadosout,X,Y,Height,Width,sField);
  fclose(dadosout);

  sprintf(filename,"%s/labels-p%3d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUlabels(dadosout,X,Y,Height,Width,label);
  fclose(dadosout);

  sprintf(filename,"%s/presence-p%3d-%d.txt",folder,index,n);
  dadosout = fopen(filename,"w");
  fprintUpresence(dadosout,X,Y,Height,Width,label);
  fclose(dadosout);
}

//----------------------------------------------------------------------

void printVorticesAndMoments(int Height,int Width, double *X,double *Y,int n,char *folder,
	                         int index,int nCnect,double *vCatalog,double *rCatalog)
{
  char filename[512+1];
  FILE *dadosout,*dadosVout;

  sprintf(filename,"%s/vortices-p%3d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");   
  err=fprintVortex(dadosVout,n,nCnect,vCatalog);
  if(err!=0){printf("problems vortices\n"); return -6;}
  fclose(dadosVout);
      
  sprintf(filename,"%s/vorticesSafe-p%3d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");
  err=fprintSafeVortex(dadosVout,n,nCnect,vCatalog,Height,Width,X,Y);
  if(err!=0){printf("problems vorticesSafe\n"); return -6;}
  fclose(dadosVout);
      
  sprintf(filename,"%s/vortexSafeMoments-p%3d-%.4f.txt",folder,index,t);
  dadosVout = fopen(filename,"w");
  err=fprintSafeVortexMoments(dadosVout,n,8,nCnect,rCatalog,Height,Width,X,Y);
  if(err!=0){printf("problems vortexSafeMoments\n"); return -6;}
  fclose(dadosVout);
}


//----------------------------------------------------------------------

  for(n=0;n<Nsnapshots;n+=1){

    t=t0+((double)n)*dt;
    if(n%10 == 0){
      printf("%d timesteps processed\n",n);
      fflush(vortexFile);
      fflush(totalVortices);
    }

    for(l=0;l<planeNum;l+=1){
      
      for(i=0;i<2*Height*Width;i+=1)
        uField[i]=0.;
      for(i=0;i<Height*Width;i+=1)
        label[i]=-1;

      dbgPrint(15,0);

      if(planeIndex>0){
      	//sprintf(filename,"%s/%g/U",foamFolder,t);
        //uFile = fopen(filename,"r");
        //if(uFile==NULL) printf("problems opening uFile - %d\n",n);

        //sprintf(filename,"%s/%g/p",foamFolder,t);
        //pFile = fopen(filename,"r");
        //if(pFile==NULL) printf("problems opening pFile - %d\n",n);

        //if(uFile == NULL || pFile == NULL)
        //  printf("Failed time = %g\n",t);
        
        safeFopen(uFile,"r","%s/%g/U",foamFolder,t);
        safeFopen(pFile,"r","%s/%g/p",foamFolder,t);

        dbgPrint(15,1);
      
        err=loadFields(Nx,Ny,Nz,uFile,pFile,node);
        if(err!=0) printf("Problems with loadFields\n");
      
        fclose(pFile); fclose(uFile);
         
        dbgPrint(15,2);
      
        err=sliceFoamField(Height,Width,pln[l],node,uField);
        if(err!=0)
          printf("Problem slicing openFoamIcoData\n");
      }
      else{
        
      }

      dbgPrint(15,3);

      err=foamScalarField(runType,Height,Width,padWidth,X,Y,Xbuff,Ybuff,
                          uField,uBuff,ux,uy,uxxx,uyyy,uxxy,
                          uxyy,gField,g2Field,v0y0,sField);
      if(err!=0){
        printf("Error in calcScalarField - %d\n",err);
        return err;
      }

      dbgPrint(15,4);

      err = floodFill(sField,Width,Height,eqClass,label);
      if(err!=0)
        printf("Problems in floodFill\n");

      err = renameLabels(Height,Width,label);
      if(err>0)
        nCnect=err;
      else
        printf("problems with renameLabels - %d\n",err);

      dbgPrint(15,4);
      
      if(n%10==0)
        printScalarFields(Height,Width,X,Y,n,folder,pln[l],sField,label);
      
      dbgPrint(15,5);

      err=extract012Momentsw2(Height,Width,nCnect,X,Y,sField,gField,label,vCatalog,vortSndMomMatrix);
      if(err!=0){
        printf("problems in extract012Momentsw2\n");
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
        rCatalog[8*i+0]=vCatalog[4*i+0];
        rCatalog[8*i+1]=vCatalog[4*i+1];
        rCatalog[8*i+2]=vCatalog[4*i+2];
        rCatalog[8*i+3]=vCatalog[4*i+3];
        rCatalog[8*i+4]=vortSndMomMatrix[4*i+0];
        rCatalog[8*i+5]=vortSndMomMatrix[4*i+1];
        rCatalog[8*i+6]=vortSndMomMatrix[4*i+2];
        rCatalog[8*i+7]=vortSndMomMatrix[4*i+3];
      }
  
      vortexAdaptiveQuickSort(rCatalog,nCnect,8,&greaterAbsCirculation);

      dbgPrint(17,0);
      
      err=histoIncVortex(nCnect,vCatalog,hG,hRc,ha,hb);
      if(err!=0){printf("problems\n"); return -5;}
    
      dbgPrint(18,0);

      /* Preparing for printing */
      if(n%10==0){
        //sprintf(filename,"%s/vortices-%.4f.txt",folder,t);
        //dadosVout = fopen(filename,"w");   
        //err=fprintVortex(dadosVout,n,nCnect,vCatalog);
        //if(err!=0){printf("problems vortices\n"); return -6;}
        //fclose(dadosVout);
      
        //sprintf(filename,"%s/vorticesSafe-%.4f.txt",folder,t);
        //dadosVout = fopen(filename,"w");
        //err=fprintSafeVortex(dadosVout,n,nCnect,vCatalog,Height,Width,X,Y);
        //if(err!=0){printf("problems vorticesSafe\n"); return -6;}
        //fclose(dadosVout);
      
        //sprintf(filename,"%s/vortexSafeMoments-%.4f.txt",folder,t);
        //dadosVout = fopen(filename,"w");
        //err=fprintSafeVortexMoments(dadosVout,n,8,nCnect,rCatalog,Height,Width,X,Y);
        //if(err!=0){printf("problems vortexSafeMoments\n"); return -6;}
        //fclose(dadosVout);

        printVorticesAndMoments(Height,Width,X,Y,n,folder,pln[l],nCnect,vCatalog,rCatalog);
        
      }

      err=fprintSafeVortexMoments(totalVortices,n,8,nCnect,rCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems vortexSafeMoments Total\n"); return -6;}
    
      err=fprintSafeVortex(vortexFile,n,nCnect,vCatalog,Height,Width,X,Y);
      if(err!=0){printf("problems in printing vortexfile\n"); return -6;}
    }
  } // End of Main loop