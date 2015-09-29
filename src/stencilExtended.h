
int XtoXbuff(int Width,float *X,float *Xbuff,int padWidth);

int uFieldTouBuff(int Height,int Width,float *X,float *Y,float *Xbuff,
                  float *Ybuff, float *uField,float *uBuff,int padWidth);

int UToGradUPadded(int Height,int Width,int type,float *x0,float *dx,
                        float *X,float *Y, float *uField, float *gField,
                        float *uBuff,float *Xbuff,float *Ybuff);