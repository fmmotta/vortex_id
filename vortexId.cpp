/*
 * File with the vortex identification prototype
 * including the floodfill algorithm applied to a
 * scalar field
 * 
 * Algorithm choice is the 2 pass flood fill algorithm
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
// #include "mt64.h"


  float sField0[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           0.,1.,1.,0.,0.,0.,0.,1.,0.,1.,
                           1.,1.,0.,1.,0.,0.,1.,1.,1.,1.,
                           0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,
                           1.,1.,1.,0.,0.,0.,0.,0.,1.,1.,
                           1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

  float sField1[] = {1.,0.,1.,0.,1.,0.,0.,0.,0.,0.,
                           1.,0.,1.,1.,0.,0.,0.,1.,0.,1.,
                           1.,0.,0.,1.,0.,0.,1.,1.,1.,1.,
                           0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
                           1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,
                           1.,0.,1.,0.,0.,1.,0.,0.,1.,1.,
                           1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           0.,0.,1.,0.,0.,0.,1.,0.,0.,0. };

  float sField2[] = {1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                           1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                           1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                           1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,
                           1.,0.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,1.,0.,1.,0.,0.,0. };
  
  float sField3[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                           1.,0.,1.,1.,0.,1.,0.,0.,0.,1.,
                           1.,0.,0.,1.,0.,0.,0.,0.,1.,1.,
                           0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                           0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
                           1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,
                           1.,0.,1.,0.,0.,1.,0.,0.,1.,1.,
                           1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           0.,0.,1.,0.,0.,0.,1.,0.,0.,0. };

  float sField4[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,1.,
                           1.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                           1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                           0.,1.,0.,1.,0.,1.,0.,0.,1.,1.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           1.,1.,1.,1.,1.,0.,0.,1.,1.,1.,
                           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           0.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                           0.,1.,0.,1.,0.,1.,0.,0.,0.,1.,
                           0.,1.,0.,1.,0.,1.,0.,0.,0.,1. };

  float sField5[] = {1.,0.,1.,0.,1.,0.,1.,0.,0.,1.,
                           1.,0.,1.,0.,1.,0.,1.,0.,1.,1.,
                           1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,
                           0.,1.,0.,0.,0.,0.,1.,0.,1.,1.,
                           0.,0.,1.,0.,0.,1.,0.,0.,0.,0.,
                           1.,0.,1.,1.,1.,0.,0.,1.,1.,1.,
                           1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                           1.,0.,1.,0.,0.,0.,1.,0.,1.,1.,
                           0.,1.,0.,1.,0.,1.,0.,0.,0.,1.,
                           0.,1.,0.,1.,0.,1.,0.,0.,0.,1. };

using namespace std;

#define _USE_MATH_DEFINES
#define dbg_prints(label,i) printf(label" debug %d\n",(i));
#define dbg_printp(label,i,j) printf(label " debug (%d,%d)\n",(i),(j));

#define bound_check(a,b,w,h) ((a)>=0)&&((a)<(h))&&((b)>=0)&&((b)<(w))

int check_neighbours(int i,int j,int *label,int Width,int Height,
                    int *nbList){
  int neighbours=0;

  // 4-way conectivity

  if(bound_check(i+1,j,Width,Height) && (label[(i+1)*Width+j]>=0)){
    nbList[neighbours*2+0] = i+1;
    nbList[neighbours*2+1] = j;
    neighbours+=1;
  }

  if(bound_check(i,j+1,Width,Height) && (label[(i)*Width+(j+1)]>=0)){
    nbList[neighbours*2+0] = i;
    nbList[neighbours*2+1] = j+1;
    neighbours+=1;
  }

  if(bound_check(i-1,j,Width,Height) && (label[(i-1)*Width+(j)]>=0)){
    nbList[neighbours*2+0] = i-1;
    nbList[neighbours*2+1] = j;
    neighbours+=1;
  }

  if(bound_check(i,j-1,Width,Height) && (label[(i)*Width+(j-1)]>=0)){
    nbList[neighbours*2+0] = i;
    nbList[neighbours*2+1] = j-1;
    neighbours+=1;
  }

  // 8-way conectivity

  if(bound_check(i+1,j+1,Width,Height) && (label[(i+1)*Width+(j+1)]>=0)){
    nbList[neighbours*2+0] = i+1;
    nbList[neighbours*2+1] = j+1;
    neighbours+=1;
  }

  if(bound_check(i+1,j-1,Width,Height) && (label[(i+1)*Width+(j-1)]>=0)){
    nbList[neighbours*2+0] = i+1;
    nbList[neighbours*2+1] = j-1;
    neighbours+=1;
  }

  if(bound_check(i-1,j+1,Width,Height) && (label[(i-1)*Width+(j+1)]>=0)){
    nbList[neighbours*2+0] = i-1;
    nbList[neighbours*2+1] = j+1;
    neighbours+=1;
  }

  if(bound_check(i-1,j-1,Width,Height) && (label[(i-1)*Width+(j-1)]>=0)){
    nbList[neighbours*2+0] = i-1;
    nbList[neighbours*2+1] = j-1;
    neighbours+=1;
  }

  return neighbours;
}

int findEq(int element,int *eqList,int pop){
  int i;
  for(i=0;i<pop;i+=1)
    if(eqList[i]==element)
      return i;
  return -1;
}

#define NumCls 1024

int checkEqClass(int eqClass[][NumCls],int eqPop[],int counter){
  int i,j,k,found;
  
  for(i=0;i<counter;i+=1){
    for(j=0;j<counter;j+=1){
      for(k=0;k<counter;k+=1){
        found = min(findEq(j,eqClass[i],eqPop[i]), findEq(k,eqClass[j],eqPop[j]));
        if(found<0)
          return 1;
      }
    }
  }
  
  return 0;
}

int floodFill(float *sField,int Width,int Height,int *label){
  int i,j,k,kpop,klabel,counter=0;
  int found,err,neighbours,nbList[8],minLabel,label2k;
  int eqClass[NumCls][NumCls]; // to be revisited
  int eqPop[NumCls];

  for(i=0;i<NumCls;i+=1){
    eqPop[i]= 0;
    for(j=0;j<NumCls;j+=1)
      eqClass[i][j] = 0;
  }

  for(i=0;i<Height*Width;i+=1)
    label[i]=-1;

  //Main flood fill loop - 1st pass
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1){

      if(sField[i*Width+j]>0){
        neighbours = check_neighbours(i,j,label,Width,Height,nbList);
        if(neighbours>0){ // need to deal with neighbours
          minLabel=label[nbList[0]*Width+nbList[1]];
          for(k=1;k<neighbours;k+=1)
            minLabel = min(label[nbList[2*k+0]*Width+nbList[2*k+1]],minLabel);
          
          label[i*Width+j]=minLabel;          
          for(k=0;k<neighbours;k+=1){           
            // check if k-th label is equivalent to minLabel
            label2k = label[nbList[2*k+0]*Width+nbList[2*k+1]];
            found = findEq(label2k,eqClass[minLabel],eqPop[minLabel]);

            // if k-th not previously found to be equivalent to minLabel
            if(found<0){
              // set k-th neighbour equivalent to minLabel  
              eqPop[minLabel]+=1;
              eqClass[minLabel][eqPop[minLabel]-1]=label2k;
            }

            // check if minLabel is equivalent to k-th
            found=findEq(minLabel,eqClass[label2k],eqPop[label2k]);
            // if minLabel is not yet equivalent to k-th neighbour
            if(found<0){
              eqPop[label2k] += 1;
              eqClass[label2k][eqPop[label2k]-1]=minLabel;
            }
          }
        }
        else{ // no neighbours - simply set new label
          label[i*Width+j]=counter;
          eqPop[counter] += 1; /* Set new label equivalent to itself */
          eqClass[counter][0]=counter;
          counter+=1;
        }
      }
    }
  }

  /* Necessary Intermediary step to assure 
   * that we have equivalence classes 
   * which satisfy transitivity
   * 
   * Though not necessary, I'm iterating 
   * because i believe that is conceptually 
   * necessary. Test necessity during tests.
   */ 
  do{
    err=0;
    for(i=0;i<counter;i+=1){
      for(j=0;j<counter;j+=1){
        for(k=0;k<counter;k+=1){
          found = min(findEq(j,eqClass[i],eqPop[i]), findEq(k,eqClass[j],
                                                            eqPop[j]));
          if(found>0){
            found = findEq(k,eqClass[i],eqPop[i]);
            if(found < 0){
              err=1;
              eqPop[i]+=1;
              eqPop[k]+=1;
              eqClass[i][eqPop[i]-1]=k;
              eqClass[k][eqPop[k]-1]=i;
            }
          }
        }
      }
    }
  }while(err==0);
  
  // this is bound to be changed
  //if(checkEqClass(eqClass,eqPop,counter) != 0) printf("problems"); 

  //Main flood fill loop - 2nd pass
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      if(label[i*Width+j]>=0){
        found = label[i*Width+j];
        minLabel=eqClass[found][0];
        for(k=1;k<eqPop[found];k+=1)
          minLabel=min(minLabel,eqClass[found][k]);
        label[i*Width+j]=minLabel;
      }

  return 0;
}

int checkFourStencil(int i,int j,int Width,int Height){
  int neighbours=0;

  // 4-way stencil check
  if(bound_check(i+1,j,Width,Height))
    neighbours+=1;
  
  if(bound_check(i,j+1,Width,Height))
    neighbours+=1;

  if(bound_check(i-1,j,Width,Height))
    neighbours+=1;
  
  if(bound_check(i,j-1,Width,Height))
    neighbours+=1;
  
  return neighbours;
}

#define ux(i,j,w) (((i)*(w)+(j))*2+0)
#define uy(i,j,w) (((i)*(w)+(j))*2+1)

int updateFrom2DVelocityField(int Width,int Height,
                              float dx,float dy,
                              float *uField,float *lambField){
  int neighbours,i,j;
  float gradU[4],detU,trU,lamb2;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      neighbours = checkFourStencil(i,j,Width,Height);
      if(neighbours!=4)
        continue;
      
      // using centered differences to calculate the gradient
      gradU[0*2+0] = (uField[ux(i+1,j,Width)]-uField[ux(i-1,j,Width)])/(2.*dx);
      gradU[0*2+1] = (uField[uy(i+1,j,Width)]-uField[uy(i-1,j,Width)])/(2.*dx);
      gradU[1*2+0] = (uField[ux(i,j+1,Width)]-uField[ux(i,j-1,Width)])/(2.*dy);
      gradU[1*2+1] = (uField[uy(i,j+1,Width)]-uField[uy(i,j+1,Width)])/(2.*dy);
    
      detU = gradU[0]*gradU[3]-gradU[1]*gradU[2];
      trU = gradU[0]+gradU[3];

      // \Delta = (trU)^2 - 4.*detU ; \Delta<0 ==> Imaginary eigenvalues
      lamb2 = detU-(trU*trU)/4.;
      if(lamb2>=0)
        lambField[i*Width+j] = sqrt(lamb2);
      else
        lambField[i*Width+j] = 0.; // I might have to change to -1
    }

  return 0;
}

int initZero(int Height,int Width, float **gFieldOut){
  int i,j,k;
  float *gField;

  if((*gFieldOut) == NULL){
    gField = (float*)malloc(4*Height*Width*sizeof(float));
    if(gField == NULL)
      return 1;
  }
  else
    gField = *gFieldOut;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gField[4*(i*Width+j)+0] = 0.;
      gField[4*(i*Width+j)+1] = 0.;
      gField[4*(i*Width+j)+2] = 0.;
      gField[4*(i*Width+j)+3] = 0.;
    }

  *gFieldOut = gField;
  return 0;
}

int addSingleOseen(int nVortex,float *parVortex, float *x0, float *dx, 
                   int Height,int Width, float **gFieldOut){
  int i,j,k;
  float gradU[2][2],*gField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;

  if(*gFieldOut==NULL)
    return 1;
  gField = *gFieldOut;
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{          
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      gField[4*(i*Width+j)+0]=gradU[0][0];
      gField[4*(i*Width+j)+1]=gradU[0][1];
      gField[4*(i*Width+j)+2]=gradU[1][0];
      gField[4*(i*Width+j)+3]=gradU[1][1];
    }

  return 0;
}

int initLambOseen2D(int nVortex,float *parVortex,
                    float *x0, float *dx, int Height,int Width,
                    float **sFieldOut){
  int i,j,k;
  float gradU[2][2],*sField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  
  if((*sFieldOut) == NULL){
    sField = (float*)malloc(Height*Width*sizeof(float));
    if(sField == NULL)
      return 1;
  }
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{          
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (gradU[0][0]*gradU[1][1]-gradU[0][1]*gradU[1][0]);
      lamb-= ((gradU[0][0]+gradU[1][1])*(gradU[0][0]+gradU[1][1]))/4.;
      
      if(lamb>0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }

  *sFieldOut = sField;
  return 0;
}

int initOseenShear2D(int nVortex,float *parVortex,
                     float *x0, float *dx, int Height,int Width,
                     float v0y0, float **sFieldOut){
  int i,j,k;
  float gradU[2][2],*sField;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  
  if((*sFieldOut) == NULL){
    sField = (float*)malloc(Height*Width*sizeof(float));
    if(sField == NULL)
      return 1;
  }
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gradU[0][1] = gradU[1][0] = gradU[1][1] = 0.;
      
      x = x0[0] + i*dx[0];
      y = x0[1] + j*dx[1];
      for(k=0;k<nVortex;k+=1){
        G = parVortex[4*k+0]; R = parVortex[4*k+1];
        a = parVortex[4*k+2]; b = parVortex[4*k+3];
        
        r2 = (x-a)*(x-a)+(y-b)*(y-b);
        r = sqrt(r2);
        
        // maybe add an if clause for small r2/(R*R)
        if(r<=0){
          gradU[0][0]=0.;
          gradU[1][1]=0.;
          gradU[0][1]=-G/(2.*M_PI*R*R);
          gradU[1][0]=G/(2.*M_PI*R*R);
        }
        else if((r>0)&&(r/R<cutoff)){
          fa=(G/(M_PI*R*R))*(0.5-0.25*(r/R)*(r/R));
          fb=(G/(M_PI*R*R))*(0.5-0.75*(r/R)*(r/R));
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
        else{
          fa = (G/(2.*M_PI*r2))*(1. - exp(-r2/(R*R)));
          fb = -1. + (1.+r2/(R*R))*exp(-r2/(R*R));
          fb *= G/(M_PI*r2);
        
          gradU[0][0] += -((x-a)*(y-b)*fb)/r2; 
          gradU[0][1] += -fa-((y-b)*(y-b)*fb)/r2;
          gradU[1][0] +=  fa+((x-a)*(x-a)*fb)/r2;
          gradU[1][1] +=  ((x-a)*(y-b)*fb)/r2;
        }
      }
      
      gradU[0][1] += v0y0; // got to improove this later
      
      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (gradU[0][0]*gradU[1][1]-gradU[0][1]*gradU[1][0]);
      lamb-= ((gradU[0][0]+gradU[1][1])*(gradU[0][0]+gradU[1][1]))/4.;
      
      if(lamb>0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }
    
  *sFieldOut = sField;
  return 0;
}

int main(int argc,char **argv){
  const int Width = 10, Height = 10, Pop=10,nVortex=3;
  int i,j,err,ngbr,found;
  int nbList[8],label[Width*Height],eqList[Pop];
  float parVortex[4*nVortex],x0[2],dx[2],xf[2],*sField=NULL;
  float x,y,v0y0 = 0.05;
                                
  err = floodFill(sField5,Width,Height,label);

  printf("\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2.0lf ",sField5[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%2d ",label[i*Width+j]+1);
    printf("\n");
  }
  

  //x0[0]=-4.75; xf[0]= 5.25; dx[0] = (xf[0]-x0[0])/Height;
  //x0[1]=-4.75; xf[1]= 5.25; dx[1] = (xf[1]-x0[1])/Width;
  
  /*
  x0[0]=-5.; xf[0]= 5.; dx[0] = (xf[0]-x0[0])/Height;
  x0[1]=-5.; xf[1]= 5.; dx[1] = (xf[1]-x0[1])/Width;
  parVortex[0]=1.; parVortex[1]=1.; parVortex[2]=-2.; parVortex[3]=0.;
  parVortex[4+0]=1.; parVortex[4+1]=1.; parVortex[4+2]=2.; parVortex[4+3]=0.;
  parVortex[8+0]=1.; parVortex[8+1]=1.; parVortex[8+2]=0.; parVortex[8+3]=4.;
  err = initOseenShear2D(nVortex,parVortex,x0,dx,Height,Width,v0y0,&sField);

  {
    FILE *dadosout;
    dadosout=fopen("initLambOseen2D.txt","w");
    for(i=0;i<Height;i+=1)
      for(j=0;j<Width;j+=1){
        x = x0[0] + i*dx[0];
        y = x0[1] + j*dx[1];
        
        fprintf(dadosout,"%f %f %f\n",x,y,sField[i*Width+j]);
      }

    fclose(dadosout);dadosout=NULL;
  }*/
  
  if(sField!=NULL)
    free(sField);
  return 0;
}
