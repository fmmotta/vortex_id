#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int XtoXbuffMirror(int Width,double *X,double *Xbuff,int padWidth){
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
    Xbuff[i] = X[0]-(X[padWidth-i]-X[0]);

  for(i=padWidth;i<Width+padWidth;i+=1)
    Xbuff[i] = X[i-padWidth];

  for(i=Width+padWidth;i<Width+2*padWidth;i+=1)
    Xbuff[i] = X[Width-1]+(X[Width-1]-X[2*(Width-1)+padWidth-i]);
 
  return 0;
}

int XtoXbuffPeriodic(int Width,double *X,double *Xbuff,int padWidth){
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
    Xbuff[i] = X[0]-(X[Width]-X[Width-i]);

  for(i=padWidth;i<Width+padWidth;i+=1)
    Xbuff[i] = X[i-padWidth];

  for(i=Width+padWidth;i<Width+2*padWidth;i+=1)
    Xbuff[i] = X[Width-1]+(X[i-(Width+padWidth)]-X[0]);
 
  return 0;
}

/*
 * Field to uBuff using mirror buffering, which changes the
 * relative weight around the borders.
 * therefore, the relative weight is something like the following:
 *  1 2 2 2 
 *  2 3 3 3
 *  2 3 3 3
 *  2 3 3 3
 */

int uFieldTouBuffMirror(int Height,int Width,double *uField,
                        double *uBuff,int padWidth){
  const int wdx=Width+2*padWidth,pdw=padWidth;
  int i,j;
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -3;
  if(uField==NULL)
    return -4;
  if(padWidth+1 >= (Width/2)-2)
    return -4;
  
  for(i=0;i<padWidth;i+=1){
    // Region 1 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(pdw-j))+1];
    }
    // Region 2 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(j-pdw))+1];
    }    
    // Region 3
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(2*(Width-1)+pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(2*(Width-1)+pdw-j))+1];
    }
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    // Region 4
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(pdw-j))+1];
    }
    // Region 5 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(j-pdw))+1];
    }    
    // Region 6 
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+1];
    }
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    // Region 7 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(pdw-j))+1];
    }
    // Region 8 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-pdw))+1];
    }    
    // Region 9
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(2*(Width-1)+pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(2*(Width-1)+pdw-j))+1];
    }
  }

  return 0;
}

int uFieldTouBuffMirrorXperiodic(int Height,int Width,double *uField,
                                 double *uBuff,int padWidth){
  const int wdx=Width+2*padWidth,pdw=padWidth;
  int i,j;
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -3;
  if(uField==NULL)
    return -4;
  if(padWidth+1 >= (Width/2)-2)
    return -4;

  for(i=0;i<padWidth;i+=1){
    // Region 1 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(Width-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(Width-j))+1];
    }
    // Region 2 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(j-pdw))+1];
    }    
    // Region 3
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((pdw-i)*Width+(j-Width-padWidth))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((pdw-i)*Width+(j-Width-padWidth))+1];
    }
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    // Region 4
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(Width-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(Width-j))+1];
    }
    // Region 5 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(j-pdw))+1];
    }    
    // Region 6 
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+1];
    }
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    // Region 7 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(Width-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(Width-j))+1];
    }
    // Region 8 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-pdw))+1];
    }    
    // Region 9
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-Width-padWidth))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((2*(Width-1)+pdw-i)*Width+(j-Width-padWidth))+1];
    }
  }

  return 0;
}

// ----------------------------------------------------------------------
//
//

int uFieldTouBuffYzeroXperiodic(int Height,int Width,double *uField,
                                double *uBuff,int padWidth){
  const int wdx=Width+2*padWidth,pdw=padWidth;
  int i,j;
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -3;
  if(uField==NULL)
    return -4;
  if(padWidth+1 >= (Width/2)-2)
    return -4;

  for(i=0;i<padWidth;i+=1){
    // Region 1 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
    // Region 2 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }    
    // Region 3
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    // Region 4
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(Width-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(Width-j))+1];
    }
    // Region 5 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(j-pdw))+1];
    }    
    // Region 6 
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(2*(Width-1)+pdw-j))+1];
    }
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    // Region 7 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
    // Region 8 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }    
    // Region 9
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
  }

  return 0;
}

// ----------------------------------------------------------------------
//
//

int uFieldTouBuffYzeroXzero(int Height,int Width,double *uField,
                                double *uBuff,int padWidth){
  const int wdx=Width+2*padWidth,pdw=padWidth;
  int i,j;
  
  if(Width<0 || Height<0)
    return -1;
  if(uBuff==NULL)
    return -3;
  if(uField==NULL)
    return -4;
  if(padWidth+1 >= (Width/2)-2)
    return -4;

  for(i=0;i<padWidth;i+=1){
    // Region 1 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
    // Region 2 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }    
    // Region 3
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
  }
  for(i=padWidth;i<Height+padWidth;i+=1){
    // Region 4
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
    // Region 5 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = uField[2*((i-pdw)*Width+(j-pdw))+0];
      uBuff[2*(i*wdx+j)+1] = uField[2*((i-pdw)*Width+(j-pdw))+1];
    }    
    // Region 6 
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
  }
  for(i=Height+padWidth;i<Height+2*padWidth;i+=1){
    // Region 7 
    for(j=0;j<padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
    // Region 8 
    for(j=padWidth;j<Width+padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }    
    // Region 9
    for(j=Width+padWidth;j<Width+2*padWidth;j+=1){
      uBuff[2*(i*wdx+j)+0] = 0.0;
      uBuff[2*(i*wdx+j)+1] = 0.0;
    }
  }

  return 0;
}


// just your run to the mill mask convolution.

int gaussianFilterUniform(int Height,int Width,int padWidth,const double *mask,
                          double *uBuff,double *uFilt)
{
  int i,j,ik,jk,idx,jdx,wdx=Width+2*padWidth,mwd=2*padWidth+1;
  
  if(Width<0 || Height<0)
    return -1;
  if(padWidth<=0)
    return -2;
  if(uBuff==NULL)
    return -3;
  if(uFilt==NULL)
    return -4;
  if(mask==NULL)
    return -5;
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uFilt[2*(i*Width+j)+0]=0.;
      uFilt[2*(i*Width+j)+1]=0.;
      idx = i+padWidth;
      jdx = j+padWidth;
      for(ik=-padWidth;ik<=padWidth;ik+=1)
        for(jk=-padWidth;jk<=padWidth;jk+=1){
          uFilt[2*(i*Width+j)+0] += mask[(ik+padWidth)*mwd+(jk+padWidth)]
                                  *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];
          uFilt[2*(i*Width+j)+1] += mask[(ik+padWidth)*mwd+(jk+padWidth)]
                                  *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }
    }
  
  return 0;
}

// Important : buffer width is padwidth+1, not the usual padWidth

int gaussianFilterNonUniform(int Height,int Width,int padWidth,
                             double *Xbuff,double *Ybuff,double *uBuff,
                             double sigma2,double *uFilt)
{
  // ipadWidth = integration pad width
  const int wdx=Width+2*padWidth,ipadWidth=padWidth-1;
  int i,j,ik,jk,idx,jdx,ii,jj;
  double wgt=0.,dist,x,y,u,v,norm,nm;
  
  if(Width<0 || Height<0)
    return -1;
  if(padWidth<=0)
    return -2;
  if(uBuff==NULL)
    return -3;
  if(uFilt==NULL)
    return -4;
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;

  //----------------------------------------------------------
  //
  //      i+2,j-2   i-1,j-1    i+2,j       i+2,j+1     i+2,j+2
  //     _____________________________________________ 
  //    |         |          |           |           |
  //    |         |          |           |           |
  //    | i+1,j-2 | i+1,j-1  | i+1,j     | i+1,j+1   | i+1,j+2
  //    |_________|__________|___________|___________|
  //    |         |          |           |           |
  //    |         |          |           |           |
  //    | i  ,j-2 | i  ,j-1  | i,j       | i  ,j+1   | i  ,j+2
  //    |_________|__________|___________|___________|
  //    |         |          |           |           |
  //    |         |          |           |           |
  //    | i-1,j-2 | i-1,j-1  | i-1,j     | i-1,j+1   | i-1,j+2
  //    |_________|__________|___________|___________|
  //    |         |          |           |           |
  //    |         |          |           |           |
  //    | i-2,j-2 | i-2,j-1  | i-2,j     | i-2,j+1   | i-2,j+2  
  //    |_________|__________|___________|___________|
  //       
  //        
  //----------------------------------------------------------

  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      norm = 0.;
      uFilt[2*(i*Width+j)+0]=0.;
      uFilt[2*(i*Width+j)+1]=0.;
      idx = i+padWidth;
      jdx = j+padWidth;
      x=Xbuff[jdx];
      y=Ybuff[idx];
      
      for(ik=-ipadWidth;ik<=ipadWidth;ik+=1)
        for(jk=-ipadWidth;jk<=ipadWidth;jk+=1){
          ii = idx+ik; jj = jdx+jk;
          
          // -------------------------------
          //  X+, Y+ region
          //
          //     ii+1,jj      ii+1,jj+1
          //     ___________  
          //    |           |           
          //    |           |           
          //    | ii,jj     | ii ,jj+1  
          //    |___________|
          //    *
          // -------------------------------

          wgt = (Ybuff[ii+1]-Ybuff[ii+0])*(Xbuff[jj+1]-Xbuff[jj+0])/64.;
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+0]*exp(-dist);
          v  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+1]*exp(-dist);
          nm = 9.*exp(-dist);

          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+1])*(x-Xbuff[jj+1]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+0)*wdx+(jj+1))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+0)*wdx+(jj+1))+1]*exp(-dist);
          nm+= 3.*exp(-dist);

          dist  = ((y-Ybuff[ii+1])*(y-Ybuff[ii+1])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+1)*wdx+(jj+0))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+1)*wdx+(jj+0))+1]*exp(-dist);
          nm+= 3.*exp(-dist);

          dist  = ((y-Ybuff[ii+1])*(y-Ybuff[ii+1])
                  +(x-Xbuff[jj+1])*(x-Xbuff[jj+1]))/(2.*sigma2);
          u += 1.*uBuff[2*((ii+1)*wdx+(jj+1))+0]*exp(-dist);
          v += 1.*uBuff[2*((ii+1)*wdx+(jj+1))+1]*exp(-dist);
          nm+= 1.*exp(-dist);

          u *= wgt;
          v *= wgt;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;

          
          // -------------------------------
          //  X-, Y+ region
          //
          //     ii+1,jj-1    ii+1,jj
          //     ___________  
          //    |           |           
          //    |           |           
          //    | ii,jj-1   | ii ,jj  
          //    |___________|
          //                 *
          // -------------------------------

          wgt = (Ybuff[ii+1]-Ybuff[ii+0])*(Xbuff[jj+0]-Xbuff[jj-1])/64.;
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+0]*exp(-dist);
          v  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+1]*exp(-dist);
          nm = 9.*exp(-dist);
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj-1])*(x-Xbuff[jj-1]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+0)*wdx+(jj-1))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+0)*wdx+(jj-1))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii+1])*(y-Ybuff[ii+1])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+1)*wdx+(jj+0))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+1)*wdx+(jj+0))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii+1])*(y-Ybuff[ii+1])
                  +(x-Xbuff[jj-1])*(x-Xbuff[jj-1]))/(2.*sigma2);
          u += 1.*uBuff[2*((ii+1)*wdx+(jj-1))+0]*exp(-dist);
          v += 1.*uBuff[2*((ii+1)*wdx+(jj-1))+1]*exp(-dist);
          nm+= 1.*exp(-dist);

          u *= wgt;
          v *= wgt;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          
          // -------------------------------
          //  X+, Y- region
          //
          //     ii,jj      ii,jj+1
          //    *___________  
          //    |           |           
          //    |           |           
          //    | ii-1,jj   | ii-1,jj+1  
          //    |___________|
          //
          // -------------------------------

          wgt = (Ybuff[ii+0]-Ybuff[ii-1])*(Xbuff[jj+1]-Xbuff[jj+0])/64.;
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+0]*exp(-dist);
          v  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+1]*exp(-dist);
          nm = 9.*exp(-dist);
          
          dist  = ((y-Ybuff[ii-1])*(y-Ybuff[ii-1])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii-1)*wdx+(jj+0))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii-1)*wdx+(jj+0))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+1])*(x-Xbuff[jj+1]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+0)*wdx+(jj+1))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+0)*wdx+(jj+1))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii-1])*(y-Ybuff[ii-1])
                  +(x-Xbuff[jj+1])*(x-Xbuff[jj+1]))/(2.*sigma2);
          u += 1.*uBuff[2*((ii-1)*wdx+(jj+1))+0]*exp(-dist);
          v += 1.*uBuff[2*((ii-1)*wdx+(jj+1))+1]*exp(-dist);
          nm+= 1.*exp(-dist);

          u *= wgt;
          v *= wgt;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          
          // -------------------------------
          //  X- Y- region
          //
          //     ii ,jj-1    ii ,jj
          //     ___________* 
          //    |           |           
          //    |           |           
          //    | ii-1,jj-1 | ii-1,jj  
          //    |___________|
          //     
          // -------------------------------

          wgt = (Ybuff[ii+0]-Ybuff[ii-1])*(Xbuff[jj+0]-Xbuff[jj-1])/64.;
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+0]*exp(-dist);
          v  = 9.*uBuff[2*((ii+0)*wdx+(jj+0))+1]*exp(-dist);
          nm = 9.*exp(-dist);
          
          dist  = ((y-Ybuff[ii-1])*(y-Ybuff[ii-1])
                  +(x-Xbuff[jj+0])*(x-Xbuff[jj+0]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii-1)*wdx+(jj+0))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii-1)*wdx+(jj+0))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii+0])*(y-Ybuff[ii+0])
                  +(x-Xbuff[jj-1])*(x-Xbuff[jj-1]))/(2.*sigma2);
          u += 3.*uBuff[2*((ii+0)*wdx+(jj-1))+0]*exp(-dist);
          v += 3.*uBuff[2*((ii+0)*wdx+(jj-1))+1]*exp(-dist);
          nm+= 3.*exp(-dist);
          
          dist  = ((y-Ybuff[ii-1])*(y-Ybuff[ii-1])
                  +(x-Xbuff[jj-1])*(x-Xbuff[jj-1]))/(2.*sigma2);
          u += 1.*uBuff[2*((ii-1)*wdx+(jj-1))+0]*exp(-dist);
          v += 1.*uBuff[2*((ii-1)*wdx+(jj-1))+1]*exp(-dist);
          nm+= 1.*exp(-dist);

          u *= wgt;
          v *= wgt;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          
          //------------------------------------------------------
          // Finished Calculating uFilt and normalization additive part
        }

      uFilt[2*(i*Width+j)+0] /= norm;
      uFilt[2*(i*Width+j)+1] /= norm;
    }

  return 0;
}