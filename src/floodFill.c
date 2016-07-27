#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "floodFill.h"

int fmind(int a,int b){
  if(a>=b)
    return b;
  else
    return a;
}

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

int checkEqClass(int eqClass[][NumCls],int eqPop[],int counter){
  int i,j,k,found;
  
  for(i=0;i<counter;i+=1){
    for(j=0;j<counter;j+=1){
      for(k=0;k<counter;k+=1){
        found = fmind(findEq(j,eqClass[i],eqPop[i]), findEq(k,eqClass[j],eqPop[j]));
        if(found<0)
          return 1;
      }
    }
  }
  
  return 0;
}

int floodFill(double *sField,int Width,int Height,int **eqClass,int *label){
  int i,j,k,counter=0;
  int found,err,neighbours,nbList[8],minLabel,label2k;
  //int eqClass[NumCls][NumCls]; // valgrind may be complaining about this
  //int **eqClass; // vai tomar no cú, isso deu certo ¬¬
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
            minLabel = fmind(label[nbList[2*k+0]*Width+nbList[2*k+1]],minLabel);
          
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
          if(counter>(NumCls-1))
            printf("moments of dispair -- floodFill\n");
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
   *
   * Consumes a lot of time, this sole loop
   * this can be problematic --
   * forget it, probably artefact of a bad
   * initial state
   */ 
   
  do{
    err=0;
    for(i=0;i<counter;i+=1){
      for(j=0;j<counter;j+=1){
        for(k=0;k<counter;k+=1){
          found = fmind(findEq(j,eqClass[i],eqPop[i]), findEq(k,eqClass[j],
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
  }while(err!=0);
  
  //Main flood fill loop - 2nd pass
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      if(label[i*Width+j]>=0){
        found = label[i*Width+j];
        minLabel=eqClass[found][0];
        for(k=1;k<eqPop[found];k+=1)
          minLabel=fmind(minLabel,eqClass[found][k]);
        label[i*Width+j]=minLabel;
      }

  return 0;
}

int renameLabels(int Height,int Width,int *label){
  int i,j,counter=0;
  int labelTag[NumCls];
  
  if((Height<=0)||(Width<=0))
    return -1;

  if(label==NULL)
    return -2;
  
  for(i=0;i<NumCls;i+=1)
    labelTag[i]=-1;

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      if(label[i*Width+j]>=0)
        if(labelTag[label[i*Width+j]]<0){
          labelTag[label[i*Width+j]]=counter;
          counter+=1;
        }

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1)
      if(label[i*Width+j]>=0)
        label[i*Width+j] = labelTag[label[i*Width+j]];

  return counter;
}