#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "stencilExtended.h"

// Y is between 0 and 1 (half channel)
// It's allways periodic on X (can be stream or span-wise)
//  
//     |
//     |
//     |
//     |
//  /\ |
//  || |
//  || | 
//  i  |
//     |
//     |_ _ _ _ _ _ _ _ _ _ _ _ _ _ 
//       
//        j ==> 
//

/*

Index manipulation example for the extended stencil.
Based calculation on the function XtoXbuff

Width = 200
padWidth = 2

i=0
padWidth+1+i = 3 + 0 = 3

i=1
padWidth+1+i = 3 + 1 = 4

3 4 0 1 2 3 4 5 6 7 8 ...

0: 3 4 0 1 2 
1: 4 0 1 2 3
2: 0 1 2 3 4
3: 1 2 3 4 5
...

i=202
i-3*padWidth-1 = 202 - 6 - 1 = 195

i=203
i-3*padWidth-1 = 203 - 6 - 1 = 196

193 194 195 196 197 198 199 195 196

...
195: 193 194 195 196 197
196: 194 195 196 197 198
197: 195 196 197 198 199
198: 196 197 198 199 195
199: 197 198 199 195 196

---------------------------------------------------

Width = 192
padWidth = 3

i=0
padWidth+1+i = 4 + 0 = 4

i=1
padWidth+1+i = 4 + 1 = 5

i=2
padWidth+1+i = 4 + 2 = 6

4 5 6 0 1 2 3 4 5 6 7 8 ...

0: 4 5 6 0 1 2 3  
1: 5 6 0 1 2 3 4
2: 6 0 1 2 3 4 5
3: 0 1 2 3 4 5 6
4: 1 2 3 4 5 6 7
...

i=192+3 = 195
i-3*padWidth-1 = 195 - 9 - 1 = 185

i=192+3+1 = 196
i-3*padWidth-1 = 196 - 9 - 1 = 186

i=192+3+1 = 197
i-3*padWidth-1 = 197 - 9 - 1 = 187

183 184 185 186 187 188 189 190 191 185 186 187

...
186: 183 184 185 186 187 188 189
187: 184 185 186 187 188 189 190
188: 185 186 187 188 189 190 191 
189: 186 187 188 189 190 191 185
190: 187 188 189 190 191 185 186
191: 188 189 190 191 185 186 187

*/

int XtoXbuff(int Width,float *X,float *Xbuff,int padWidth){
  int i;

  if(Width<=0)
    return 1;
  if(X==NULL)
    return 2;
  if(Xbuff==NULL)
    return 3;
  if(padWidth >= (Width/2)-2)
    return 4;

  for(i=0;i<padWidth;i+=1)
    Xbuff[i] = X[2*padWidth-i];

  for(i=padWidth;i<Width+padWidth;i+=1)
    Xbuff[i] = X[i-padWidth];

  for(i=Width+padWidth;i<Width+2*padWidth;i+=1)
    Xbuff[i] = X[i-3*padWidth-1];
 
  return 0;
}

/*
  Must test this function if a lot of care
 */

int gFieldTogBuff(int Height,int Width,float *X,float *Y,float *Xbuff,
                  float *Ybuff, float *gField,float *gBuff,int padWidth)
{
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(gBuff==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  for(i=0;i<padWidth;i+=1){
    for(j=0;j<padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(2*padWidth-i)*Width+(2*padWidth-j)];

    for(j=padWidth;j<Width+padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(2*padWidth-i)*Width+j];
    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(2*padWidth-i)*Width+(j-3*padWidth-1)];
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    for(j=0;j<padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[i*Width+(2*padWidth-j)];

    for(j=padWidth;j<Width+padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[i*Width+j];
    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[i*Width+(j-3*padWidth-1)];
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    for(j=0;j<padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(i-3*padWidth-1)*Width+(2*padWidth-j)];

    for(j=padWidth;j<Width+padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(i-3*padWidth-1)*Width+j];
    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1)
      gBuff[i*(Width+2*padWidth)+j] = gField[(i-3*padWidth-1)*Width+(j-3*padWidth-1)];
  }

  return 0;
}

int UToGradUPadded(int Height,int Width,int type,float *x0,float *dx,
                        float *X,float *Y, float *uField, float *gField,
                        float *gBuff,float *Xbuff,float *Ybuff){
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(dx==NULL)
    return -2;
  if(uField==NULL)
    return -3;
  if(gField==NULL)
    return -4;
  if(X==NULL)
    return -5;
  if(Y==NULL)
    return -6;
  


  return 0;
}