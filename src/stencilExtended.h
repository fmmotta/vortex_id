
int XtoXbuff(int Width,float *X,float *Xbuff,int padWidth);

int uFieldTouBuff(int Height,int Width,float *uField,float *uBuff,int padWidth);

int UToGrad3UPadded(int Height,int Width,float *gField,float *uBuff,float *Xbuff,float *Ybuff);

int UtoUx5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff);

int UtoUy5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff);

int UtoUxx5point(int Height,int Width,float *uDel,float *uBuff,
                 float *Xbuff,float *Ybuff);

int UtoUyy5point(int Height,int Width,float *uDel,float *uBuff,
                float *Xbuff,float *Ybuff);

int UtoUxxx5point(int Height,int Width,float *uDel,float *uBuff,
                  float *Xbuff,float *Ybuff);

int UtoUyyy5point(int Height,int Width,float *uDel,float *uBuff,
                  float *Xbuff,float *Ybuff);

int gradU2UtoLambda(int Height,int Width, float *gField,float *g2Field,
                    float **sFieldOut);