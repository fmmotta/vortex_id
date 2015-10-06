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
4: 2 3 4 5 6
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
    Xbuff[i] = X[padWidth+1+i];

  for(i=padWidth;i<Width+padWidth;i+=1)
    Xbuff[i] = X[i-padWidth];

  for(i=Width+padWidth;i<Width+2*padWidth;i+=1)
    Xbuff[i] = X[i-3*padWidth-1];
 
  return 0;
}

/*
  Must test this function if a lot of care
 */

int uFieldTouBuff(int Height,int Width,float *uField,float *uBuff,int padWidth)
{
  int i,j,ip,jp;
  float um,u0,up,vm,v0,vp;
  float ax,ay,bx,by;
  float cp,c0,cm;
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -3;
  if(uField==NULL)
    return -4;
  if(padWidth >= (Width/2)-2)
    return -4;
  
  for(i=0;i<padWidth;i+=1){
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((1+i+padWidth)*Width+(1+j+padWidth))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((1+i+padWidth)*Width+(1+j+padWidth))+1];
    }
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((padWidth+1+i)*Width+(j-padWidth))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((padWidth+1+i)*Width+(j-padWidth))+1];
    }    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((padWidth+1+i)*Width+(j-3*padWidth-1))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((padWidth+1+i)*Width+(j-3*padWidth-1))+1];
    }
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-padWidth)*Width+(1+j+padWidth))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-padWidth)*Width+(1+j+padWidth))+1];
    }
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-padWidth)*Width+(j-padWidth))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-padWidth)*Width+(j-padWidth))+1];
    }    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-padWidth)*Width+(j-3*padWidth-1))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-padWidth)*Width+(j-3*padWidth-1))+1];
    }
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-3*padWidth-1)*Width+(padWidth+1+j))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-3*padWidth-1)*Width+(padWidth+1+j))+1];
    }
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-3*padWidth-1)*Width+(j-padWidth))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-3*padWidth-1)*Width+(j-padWidth))+1];
    }
    
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*(Width+2*padWidth)+j)+0] = uField[2*((i-3*padWidth-1)*Width+(j-3*padWidth-1))+0];
      uBuff[2*(i*(Width+2*padWidth)+j)+1] = uField[2*((i-3*padWidth-1)*Width+(j-3*padWidth-1))+1];
    }
  }

  return 0;
}

/*
  Block derivative series:
  Calculate in each individual funcion a single type of derivative, and
  then 
 */

/*
 * \partial_x \vec u = (\partial_x u,\partial_x v) stencil calculation
 */
int UtoUx5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float a1,a2,a3,a4,c[MaskWidth]; // x positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;

      a1 = Xbuff[j+padWidth-2] - Xbuff[j+padWidth];
      a2 = Xbuff[j+padWidth-1] - Xbuff[j+padWidth];
      a3 = Xbuff[j+padWidth+1] - Xbuff[j+padWidth];
      a4 = Xbuff[j+padWidth+2] - Xbuff[j+padWidth];
      
      c[0] = -(a2*a3*a4)/(a1*(a1-a2)*(a1-a3)*(a1-a4));
      c[1] = -(a1*a3*a4)/(a2*(a2-a1)*(a2-a3)*(a2-a4));
      c[2] = -1./a1 -1./a2 -1./a3 -1./a4;
      c[3] = -(a1*a2*a4)/(a3*(a3-a1)*(a3-a2)*(a3-a4));
      c[4] = -(a1*a2*a3)/(a4*(a4-a1)*(a4-a2)*(a4-a3));

      uDel[2*(i*Width+j)+0] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+0];
      uDel[2*(i*Width+j)+1] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+1];

      uDel[2*(i*Width+j)+0] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+0];
      uDel[2*(i*Width+j)+1] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+1];

      uDel[2*(i*Width+j)+0] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+0];
      uDel[2*(i*Width+j)+1] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+1];

      uDel[2*(i*Width+j)+0] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+0];
      uDel[2*(i*Width+j)+1] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+1];
    }
  }

  return 0;
}

/*
 * \partial_y \vec u = (\partial_y u,\partial_y v) stencil calculation
 */ 
int UtoUy5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float b1,b2,b3,b4,d[MaskWidth]; // y positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    b1 = Ybuff[i+padWidth-2] - Ybuff[i+padWidth];
    b2 = Ybuff[i+padWidth-1] - Ybuff[i+padWidth];
    b3 = Ybuff[i+padWidth+1] - Ybuff[i+padWidth];
    b4 = Ybuff[i+padWidth+2] - Ybuff[i+padWidth];

    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;

      d[0] = -(b2*b3*b4)/(b1*(b1-b2)*(b1-b3)*(b1-b4));
      d[1] = (b1*b3*b4)/(b2*(b1-b2)*(b2-b3)*(b2-b4));
      d[2] = -1./b1 -1./b2 -(b3+b4)/(b3*b4);
      d[3] = (b1*b2*b4)/(b3*(b1-b3)*(b3-b2)*(b3-b4));
      d[4] = (b1*b2*b3)/(b4*(b1-b4)*(b4-b2)*(b4-b3));

      uDel[2*(i*Width+j)+0] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+1];
    }
  }

  return 0;
}

/*
 * \partial_xx \vec u = (\partial_xx u,\partial_xx v) stencil calculation
 */
int UtoUxx5point(int Height,int Width,float *uDel,float *uBuff,
                 float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float a1,a2,a3,a4,c[MaskWidth]; // x positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;

      a1 = Xbuff[j+padWidth-2] - Xbuff[j+padWidth];
      a2 = Xbuff[j+padWidth-1] - Xbuff[j+padWidth];
      a3 = Xbuff[j+padWidth+1] - Xbuff[j+padWidth];
      a4 = Xbuff[j+padWidth+2] - Xbuff[j+padWidth];
      
      c[0] = 2.*(a3*a4+a2*(a3+a4)); c[0] /= a1*(a1-a2)*(a1-a3)*(a1-a4);
      c[1] = 2.*(a3*a4+a1*(a3+a4)); c[1] /= a2*(a2-a1)*(a2-a3)*(a2-a4);
      c[2] = 2.*(a3*a4+a2*(a3+a4)+a1*(a2+a3+a4)); c[2] /= a1*a2*a3*a4;
      c[3] = 2.*(a2*a4+a2*(a1+a4)); c[3] /= a3*(a3-a1)*(a3-a2)*(a3-a4);
      c[4] = 2.*(a2*a3+a1*(a2+a3)); c[4] /= a4*(a1-a2)*(a1-a3)*(a1-a4);

      uDel[2*(i*Width+j)+0] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+0];
      uDel[2*(i*Width+j)+1] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+1];

      uDel[2*(i*Width+j)+0] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+0];
      uDel[2*(i*Width+j)+1] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+1];

      uDel[2*(i*Width+j)+0] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+0];
      uDel[2*(i*Width+j)+1] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+1];

      uDel[2*(i*Width+j)+0] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+0];
      uDel[2*(i*Width+j)+1] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+1];
    }
  }

  return 0;
}

/*
 * \partial_yy \vec u = (\partial_yy u,\partial_yy v) stencil calculation
 */
int UtoUyy5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float b1,b2,b3,b4,d[MaskWidth]; // y positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    b1 = Ybuff[i+padWidth-2] - Ybuff[i+padWidth];
    b2 = Ybuff[i+padWidth-1] - Ybuff[i+padWidth];
    b3 = Ybuff[i+padWidth+1] - Ybuff[i+padWidth];
    b4 = Ybuff[i+padWidth+2] - Ybuff[i+padWidth];

    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;
      
      d[0] = 2.*(b3*b4+b2*(b3+b4)); d[0] /= b1*(b1-b2)*(b1-b3)*(b1-b4);
      d[1] = 2.*(b3*b4+b1*(b3+b4)); d[1] /= b2*(b2-b1)*(b2-b3)*(b2-b4);
      d[2] = 2.*(b3*b4+b2*(b3+b4)+b1*(b2+b3+b4)); d[2] /= b1*b2*b3*b4;
      d[3] = 2.*(b2*b4+b2*(b1+b4)); d[3] /= b3*(b3-b1)*(b3-b2)*(b3-b4);
      d[4] = 2.*(b2*b3+b1*(b2+b3)); d[4] /= b4*(b1-b2)*(b1-b3)*(b1-b4);

      uDel[2*(i*Width+j)+0] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+1];
    }
  }

  return 0;
}

/*
 * \partial_xxx \vec u = (\partial_xxx u,\partial_xxx v) stencil calculation
 */
int UtoUxxx5point(int Height,int Width,float *uDel,float *uBuff,
                 float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float a1,a2,a3,a4,c[MaskWidth]; // x positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;

      a1 = Xbuff[j+padWidth-2] - Xbuff[j+padWidth];
      a2 = Xbuff[j+padWidth-1] - Xbuff[j+padWidth];
      a3 = Xbuff[j+padWidth+1] - Xbuff[j+padWidth];
      a4 = Xbuff[j+padWidth+2] - Xbuff[j+padWidth];
      
      c[0] = -6.*(a2+a3+a4); c[0] /= a1*(a1-a2)*(a1-a3)*(a1-a4);
      c[1] = -6.*(a1+a3+a4); c[1] /= a2*(a2-a1)*(a2-a3)*(a2-a4);
      c[2] = -6.*(a1+a2+a3+a4); c[2] /= a1*a2*a3*a4;
      c[3] = -6.*(a1+a2+a4); c[3] /= a3*(a3-a1)*(a3-a2)*(a3-a4);
      c[4] = -6.*(a1+a2+a3); c[4] /= a4*(a4-a1)*(a4-a2)*(a4-a3);

      uDel[2*(i*Width+j)+0] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+0];
      uDel[2*(i*Width+j)+1] += c[0]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-2))+1];

      uDel[2*(i*Width+j)+0] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+0];
      uDel[2*(i*Width+j)+1] += c[1]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth-1))+1];

      uDel[2*(i*Width+j)+0] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += c[2]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+0];
      uDel[2*(i*Width+j)+1] += c[3]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+1))+1];

      uDel[2*(i*Width+j)+0] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+0];
      uDel[2*(i*Width+j)+1] += c[4]*uBuff[2*((i+padWidth)*(Width+2*padWidth)+(j+padWidth+2))+1];
    }
  }

  return 0;
}

/*
 * \partial_yyy \vec u = (\partial_yyy u,\partial_yyy v) stencil calculation
 */
int UtoUyyy5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float b1,b2,b3,b4,d[MaskWidth]; // y positions and weights
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(uDel==NULL)
    return -7;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  // and the width of uBuff is Width+2*padWidth
  for(i=0;i<Height;i+=1){
    b1 = Ybuff[i+padWidth-2] - Ybuff[i+padWidth];
    b2 = Ybuff[i+padWidth-1] - Ybuff[i+padWidth];
    b3 = Ybuff[i+padWidth+1] - Ybuff[i+padWidth];
    b4 = Ybuff[i+padWidth+2] - Ybuff[i+padWidth];

    for(j=0;j<Width;j+=1){
      uDel[2*(i*Width+j)+0]=0.;
      uDel[2*(i*Width+j)+1]=0.;
      
      d[0] = -6.*(b2+b3+b4); d[0] /= b1*(b1-b2)*(b1-b3)*(b1-b4);
      d[1] = -6.*(b1+b3+b4); d[1] /= b2*(b2-b1)*(b2-b3)*(b2-b4);
      d[2] = -6.*(b1+b2+b3+b4); d[2] /= b1*b2*b3*b4;
      d[3] = -6.*(b1+b2+b4); d[3] /= b3*(b3-b1)*(b3-b2)*(b3-b4);
      d[4] = -6.*(b1+b2+b3); d[4] /= b4*(b4-b1)*(b4-b2)*(b4-b3);

      uDel[2*(i*Width+j)+0] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[0]*uBuff[2*((i+padWidth-2)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[1]*uBuff[2*((i+padWidth-1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[2]*uBuff[2*((i+padWidth+0)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[3]*uBuff[2*((i+padWidth+1)*(Width+2*padWidth)+(j+padWidth))+1];

      uDel[2*(i*Width+j)+0] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+0];
      uDel[2*(i*Width+j)+1] += d[4]*uBuff[2*((i+padWidth+2)*(Width+2*padWidth)+(j+padWidth))+1];
    }
  }

  return 0;
}


/*
 Danger: I'm not sure if this direct product of stencils is valid!
 Obs: this function is failing miserably.
 */
int UToGrad3UPadded(int Height,int Width,float *gField,float *uBuff,
                    float *Xbuff,float *Ybuff)
{
  const int MaskWidth=5,padWidth=2;
  int i,j,ip,jp,ii,jj;
  float a1,a2,a3,a4,c[MaskWidth]; // x positions and weights
  float b1,b2,b3,b4,d[MaskWidth]; // y-positions and weights
  float uBlock[MaskWidth][MaskWidth][2], stencil[MaskWidth][MaskWidth]; // local block and stencil
  
  if(Width<0 || Height<0)
    return -1;
  if(gField==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(padWidth >= (Width/2)-2)
    return 4;
  
  // Need to add +padWidth in order to offset the 
  // value of the padding. the reference value is
  // thus i+padWidth and not i alone, the same for j
  for(i=0;i<Height;i+=1){
    b1 = Ybuff[i+padWidth-2] - Ybuff[i+padWidth];
    b2 = Ybuff[i+padWidth-1] - Ybuff[i+padWidth];
    b3 = Ybuff[i+padWidth+1] - Ybuff[i+padWidth];
    b4 = Ybuff[i+padWidth+2] - Ybuff[i+padWidth];

    for(j=0;j<Width;j+=1){
      gField[4*(i*Width+j)+0]=0.;
      gField[4*(i*Width+j)+1]=0.;
      gField[4*(i*Width+j)+2]=0.;
      gField[4*(i*Width+j)+3]=0.;

      a1 = Xbuff[j+padWidth-2] - Xbuff[j+padWidth];
      a2 = Xbuff[j+padWidth-1] - Xbuff[j+padWidth];
      a3 = Xbuff[j+padWidth+1] - Xbuff[j+padWidth];
      a4 = Xbuff[j+padWidth+2] - Xbuff[j+padWidth];

      for(ii=0;ii<=2*padWidth;ii+=1)
        for(jj=0;jj<=2*padWidth;jj+=1){
          uBlock[ii][jj][0] = uBuff[2*( (i+ii)*(Width+2*padWidth)+(j+jj) )+0];
          uBlock[ii][jj][1] = uBuff[2*( (i+ii)*(Width+2*padWidth)+(j+jj) )+1];
        }
      
      /******************************************************************/

      // \partial_x^2
      c[2]  = a3*a4 + a2*(a3+a4) + a1*(a2+a3+a4);
      c[2] *= 2./(a1*a2*a3*a4); // c0
      c[0]  = a3*a4+a2*(a3+a4);
      c[0] *= 2./(a1*(a1-a2)*(a1-a3)*(a1-a4)); // c1
      c[1]  = a3*a4+a1*(a3+a4);
      c[1] *= 2./(a2*(a2-a1)*(a2-a3)*(a2-a4)); // c2
      c[3]  = a2*a4+a1*(a2+a4);
      c[3] *= 2./(a3*(a3-a1)*(a3-a2)*(a3-a4)); // c3
      c[4]  = a2*a3+a1*(a2+a3);
      c[4] *= 2./(a4*(a4-a1)*(a4-a2)*(a4-a3)); // c4

      // \partial_y 
      d[2] = -1./b1 -1./b2 - (b3+b4)/(b3*b4); // d0
      d[0] = -b2*b3*b4/(b1*(b1-b2)*(b1-b3)*(b1-b4)); // d1
      d[1] = -b1*b3*b4/(b2*(b2-b1)*(b2-b3)*(b2-b4)); // d2
      d[3] = -b1*b2*b4/(b3*(b3-b1)*(b3-b2)*(b3-b4)); // d3
      d[4] = -b1*b2*b3/(b4*(b4-b1)*(b4-b2)*(b4-b3)); // d4

      // \partial_x^2 \partial_y u & \partial_x^2 \partial_y^2 v      
      for(ii=0;ii<MaskWidth;ii+=1)
        for(jj=0;jj<MaskWidth;jj+=1){
          gField[4*(i+Width+j)+0] +=  d[ii]*c[jj]*uBlock[ii][jj][1]; // + v_xxy
          gField[4*(i+Width+j)+3] += -d[ii]*c[jj]*uBlock[ii][jj][1]; // - v_xxy
          gField[4*(i+Width+j)+2] +=  d[ii]*c[jj]*uBlock[ii][jj][0]; // + u_xxy
        }
      
      /******************************************************************/
      
      // \partial_x 
      c[2] = -1./a1 -1./a2 - (a3+a4)/(a3*a4); // c0
      c[0] = -a2*a3*a4/(a1*(a1-a2)*(a1-a3)*(a1-a4)); // c1
      c[1] = -a1*a3*a4/(a2*(a2-a1)*(a2-a3)*(a2-a4)); // c2
      c[3] = -a1*a2*a4/(a3*(a3-a1)*(a3-a2)*(a3-a4)); // c3
      c[4] = -a1*a2*a3/(a4*(a4-a1)*(a4-a2)*(a4-a3)); // c4
      
      // \partial_y^2
      d[2]  = b3*b4 + b2*(b3+b4) + b1*(b2+b3+b4);
      d[2] *= 2./(b1*b2*b3*b4); // d0
      d[0]  = b3*b4+b2*(b3+b4);
      d[0] *= 2./(b1*(b1-b2)*(b1-b3)*(b1-b4)); // d1
      d[1]  = b3*b4+b1*(b3+b4);
      d[1] *= 2./(b2*(b2-b1)*(b2-b3)*(b2-b4)); // d2
      d[3]  = b2*b4+b1*(b2+b4);
      d[3] *= 2./(b3*(b3-b1)*(b3-b2)*(b3-b4)); // d3
      d[4]  = b2*b3+b1*(b2+b3);
      d[4] *= 2./(b4*(b4-b1)*(b4-b2)*(b4-b3)); // d4
      
      // \partial_x \partial_y^2 u & \partial_x \partial_y^2 v      
      for(ii=0;ii<MaskWidth;ii+=1)
        for(jj=0;jj<MaskWidth;jj+=1){
          gField[4*(i+Width+j)+0] += -d[ii]*c[jj]*uBlock[ii][jj][0]; // - u_xyy
          gField[4*(i+Width+j)+1] +=  d[ii]*c[jj]*uBlock[ii][jj][1]; // + v_xyy
          gField[4*(i+Width+j)+3] +=  d[ii]*c[jj]*uBlock[ii][jj][0]; // + u_xxy
        }
      
      /******************************************************************/
      
      // \partial_x^3 
      c[2] = - (6.*(a1+a2+a3+a4))/(a1*a2*a3*a4);
      c[0] = - (6.*(a2+a3+a4))/(a1*(a1-a2)*(a1-a3)*(a1-a4));
      c[1] = - (6.*(a1+a3+a4))/(a2*(a2-a1)*(a2-a3)*(a2-a4));
      c[3] = - (6.*(a1+a2+a4))/(a3*(a3-a1)*(a3-a2)*(a3-a4));
      c[4] = - (6.*(a1+a2+a3))/(a4*(a4-a1)*(a4-a2)*(a4-a3));

      // \partial_y^3 
      d[2] = - (6.*(b1+b2+b3+b4))/(b1*b2*b3*b4);
      d[0] = - (6.*(b2+b3+b4))/(b1*(b1-b2)*(b1-b3)*(b1-b4));
      d[1] = - (6.*(b1+b3+b4))/(b2*(b2-b1)*(b2-b3)*(b2-b4));
      d[3] = - (6.*(b1+b2+b4))/(b3*(b3-b1)*(b3-b2)*(b3-b4));
      d[4] = - (6.*(b1+b2+b3))/(b4*(b4-b1)*(b4-b2)*(b4-b3));

      // \partial_y^3 u
      for(ii=0;ii<MaskWidth;ii+=1)
        gField[4*(i+Width+j)+1] += - d[ii]*uBlock[ii][padWidth][0]; // -u_yyy

      // \partial_x^3 v
      for(jj=0;jj<MaskWidth;jj+=1)
        gField[4*(i+Width+j)+2] += - c[jj]*uBlock[padWidth][jj][1]; // -v_xxx
    }
  }

  return 0;
}

int gradU2UtoLambda(int Height,int Width, float *gField,float *g2Field,
                    float **sFieldOut)
{
  int i,j,k;
  float gradU[2][2],grad2U[2][2],*sField,w,D2w;
  float a,b,G,R,x,y,fa,fb,r2,r,lamb,cutoff=0.001;
  
  if((*sFieldOut) == NULL) 
    return -10;
  else
    sField = *sFieldOut;
    
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      gradU[0][0] = gField[4*(i*Width+j)+0];
      gradU[0][1] = gField[4*(i*Width+j)+1];
      gradU[1][0] = gField[4*(i*Width+j)+2];
      gradU[1][1] = gField[4*(i*Width+j)+3];

      w = gradU[0][1]-gradU[1][0];

      grad2U[0][0] = g2Field[4*(i*Width+j)+0];
      grad2U[0][1] = g2Field[4*(i*Width+j)+1];
      grad2U[1][0] = g2Field[4*(i*Width+j)+2];
      grad2U[1][1] = g2Field[4*(i*Width+j)+3];

      D2w = grad2U[0][1]-grad2U[1][0];

      // \Delta = (tr gU)^2-4.*det gU; \Delta<0 ==> Imaginary eigenvalue
      // (lamb)^2 = - 4.*\Delta;
      lamb =  (grad2U[0][0]*grad2U[1][1]-grad2U[0][1]*grad2U[1][0]);
      lamb-= ((grad2U[0][0]+grad2U[1][1])*(grad2U[0][0]+grad2U[1][1]))/4.;
      
      if(lamb>0. && (w*D2w)<0.)
        sField[i*Width+j] = sqrt(lamb);
      else
        sField[i*Width+j] = 0.;
    }
    
  *sFieldOut = sField;

  return 0;
}
