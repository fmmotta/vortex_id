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

int floodFill(float *sField,int Width,int Height,int *label){
  int i,j,k,counter=0;
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

        printf("\n");
      }
    }
  }

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

/*
int initBinaryRandomField(unsigned long long int seed,
                          float smin,float smax,
                          int Width,int Height,float **sField){
  int i,j;
  long int ch;
  float *sfield;
  sfield = (float*)malloc(Width*Height*sizeof(float));
  if(sfield==NULL)
    return 1;

  init_genrand64(seed);

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      ch = genrand64_int64();
      if(ch%2 == 0)
        sfield[i*Width+j] = smin;
      else
        sfield[i*Width+j] = smax;
    }      

  *sField = sfield;
  return 0;
}*/

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

      // \Delta = (trU)^2 - 4.*detU ; \Delta < 0 ==> Imaginary eigenvalues
      lamb2 = detU-(trU*trU)/4.;
      if(lamb2>=0)
        lambField[i*Width+j] = sqrt(lamb2);
      else
        lambField[i*Width+j] = 0.; // I might have to change to -1
    }

  return 0;
}

int main(int argc,char **argv){
  const int Width = 10, Height = 10, Pop=10;
  int i,j,err,ngbr,found;
  int nbList[8],label[Width*Height],eqList[Pop];
  float sField0[Width*Height] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                                 0.,1.,1.,0.,0.,0.,0.,1.,0.,1.,
                                 1.,1.,0.,1.,0.,0.,1.,1.,1.,1.,
                                 0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                                 1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,
                                 1.,1.,1.,0.,0.,0.,0.,0.,1.,1.,
                                 1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,
                                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

  float sField1[Width*Height] = {1.,0.,1.,0.,1.,0.,0.,0.,0.,0.,
                                 1.,0.,1.,1.,0.,0.,0.,1.,0.,1.,
                                 1.,0.,0.,1.,0.,0.,1.,1.,1.,1.,
                                 0.,1.,1.,0.,0.,0.,0.,0.,1.,0.,
                                 0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
                                 1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,
                                 1.,0.,1.,0.,0.,1.,0.,0.,1.,1.,
                                 1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,
                                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                                 0.,0.,1.,0.,0.,0.,1.,0.,0.,0. };

  

  err = floodFill(sField1,Width,Height,label);

  printf("\n\nsField:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%.0lf ",sField1[i*Width+j]);
    printf("\n");
  }
  printf("\n");

  printf("\nlabel:\n");
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1)
      printf("%d ",label[i*Width+j]+1);
    printf("\n");
  }

  return 0;
}