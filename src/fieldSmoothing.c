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

// just your run to the mill mask convolution.

int gaussianFilterUniform(int Height,int Width,int padWidth,const double *mask,
                          double *Xbuff,double *Ybuff,double *uBuff,
                          double *uFilt)
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
  if(Xbuff==NULL)
    return -5;
  if(Ybuff==NULL)
    return -6;
  if(mask==NULL)
    return -7;
  
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
  const int wdx=Width+2*padWidth;
  int i,j,ik,jk,idx,jdx;
  double wgt=0.,dist,scale=1./(2.0*M_PI*sigma2);
  
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
  
  for(i=0;i<Height;i+=1)
    for(j=0;j<Width;j+=1){
      uFilt[2*(i*Width+j)+0]=0.;
      uFilt[2*(i*Width+j)+1]=0.;
      idx = i+padWidth;
      jdx = j+padWidth;
      
      // First Block - Lower block - Regions 1,2,3
      ik=-(padWidth-1);
      {
        //--- Region 1 ----
        jk = -(padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*( 9.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-0])

                           + 1.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-0]) );
          
          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 2 ----
        for(jk=-(padWidth-1)+1;jk<=(padWidth-1)-1;jk+=1){
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*(12.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 4.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-0]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 3 ----
        jk =  (padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;
          
          wgt = (1.0/64.0)*( 9.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-0])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1])

                           + 1.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }
      }

      // Second Block - Middle block - Regions 4,5,6
      for(ik=-(padWidth-1)+1;ik<=(padWidth-1)-1;ik+=1){
        //--- Region 4 ----
        jk = -(padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*(12.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 4.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 5 ----
        for(jk=-(padWidth-1)+1;jk<=(padWidth-1)-1;jk+=1){
          dist  = (Xbuff[jdx+ik]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*(16.0*fabs((Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                     *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])) );

          //printf("%d %d %d %d: %lf %lf\n",i,j,ik,jk,dist,wgt);

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 6 ----
        jk =  (padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*(12.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 4.0*(Xbuff[jdx+jk+0]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }
      }

      // Third Block - Top block - Regions 7,8,9
      ik= (padWidth-1);
      {
        //--- Region 7 ----
        jk = -(padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*( 9.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])

                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1])

                           + 1.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-0])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 8 ----
        for(jk=-(padWidth-1)+1;jk<=(padWidth-1)-1;jk+=1){
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;

          wgt = (1.0/64.0)*(12.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])
                           + 4.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }

        //--- Region 9 ----
        jk =  (padWidth-1);
        {
          dist  = (Xbuff[jdx+jk]-Xbuff[jdx])*(Xbuff[jdx+jk]-Xbuff[jdx])
                + (Ybuff[idx+ik]-Ybuff[idx])*(Ybuff[idx+ik]-Ybuff[idx]);
          dist /= 2.0*sigma2;
          
          wgt = (1.0/64.0)*( 9.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])
                           + 3.0*(Xbuff[jdx+jk+1]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1])
                           + 3.0*(Xbuff[jdx+jk+0]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+1]-Ybuff[idx+ik-1])
                           + 1.0*(Xbuff[jdx+jk+0]-Xbuff[jdx+jk-1])
                                *(Ybuff[idx+ik+0]-Ybuff[idx+ik-1]) );

          uFilt[2*(i*Width+j)+0] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+0];

          uFilt[2*(i*Width+j)+1] += wgt*exp(-dist)
                                       *uBuff[2*((idx+ik)*wdx+(jdx+jk))+1];
        }
      }
      uFilt[2*(i*Width+j)+0] *= scale;
      uFilt[2*(i*Width+j)+1] *= scale;
    }

  return 0;
}

int gaussianFilterNonUniform2(int Height,int Width,int padWidth,
                              double *Xbuff,double *Ybuff,double *uBuff,
                              double sigma2,double *uFilt)
{
  // ipadWidth = integration pad width
  const int wdx=Width+2*padWidth,ipadWidth=padWidth-1;
  int i,j,ik,jk,idx,jdx,ii,jj;
  double scale=1./(2.*M_PI*sigma2);
  double wgt=0.,dist,x,y,u,v,norm,nm;

  printf("sigma2=%lf\n",sigma2);
  printf("scale=%lf\n",scale);
  
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
        for(jk=-ipadWidth;ik<=ipadWidth;ik+=1){
          ii = idx+ik; jj = jdx+jk;
          //--------------------------------------------------------
          // Region Y+ X+ 
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

          u *= wgt;//*scale; 
          v *= wgt;//*scale;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          //--------------------------------------------------------
          // Region Y+ X- 
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

          u *= wgt;//*scale; 
          v *= wgt;//*scale;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          //--------------------------------------------------------
          // Region Y- X+ 
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

          u *= wgt;//*scale; 
          v *= wgt;//*scale;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          //--------------------------------------------------------
          // Region Y- X- 
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

          u *= wgt;//*scale; 
          v *= wgt;//*scale;
          nm*= wgt;
         
          uFilt[2*(i*Width+j)+0] += u;
          uFilt[2*(i*Width+j)+1] += v;
          norm += nm;
          //------------------------------------------------------
          // Finished Calculating uFilt -- filtered uField
        }

      //uFilt[2*(i*Width+j)+0] *= scale;
      //uFilt[2*(i*Width+j)+1] *= scale;
      uFilt[2*(i*Width+j)+0] /= norm;
      uFilt[2*(i*Width+j)+1] /= norm;
    }

  return 0;
}