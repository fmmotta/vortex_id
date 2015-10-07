
int XtoXbuff(int Width,double *X,double *Xbuff,int padWidth);

int uFieldTouBuff(int Height,int Width,double *uField,double *uBuff,int padWidth);

int UToGrad3UPadded(int Height,int Width,double *gField,double *uBuff,double *Xbuff,double *Ybuff);

int UtoUx5point(int Height,int Width,double *uDel,double *uBuff,
                double *Xbuff,double *Ybuff);

int UtoUy5point(int Height,int Width,double *uDel,double *uBuff,
                double *Xbuff,double *Ybuff);

int UtoUxx5point(int Height,int Width,double *uDel,double *uBuff,
                 double *Xbuff,double *Ybuff);

int UtoUyy5point(int Height,int Width,double *uDel,double *uBuff,
                double *Xbuff,double *Ybuff);

int UtoUxxx5point(int Height,int Width,double *uDel,double *uBuff,
                  double *Xbuff,double *Ybuff);

int UtoUyyy5point(int Height,int Width,double *uDel,double *uBuff,
                  double *Xbuff,double *Ybuff);

int gradU2UtoLambda(int Height,int Width, double *gField,double *g2Field,
                    double **sFieldOut);