int XtoXbuffMirror(int Width,double *X,double *Xbuff,int padWidth);

int uFieldTouBuffMirror(int Height,int Width,double *uField,
                        double *uBuff,int padWidth);

int gaussianFilterUniform(int Height,int Width,int padWidth,const double *mask,
                          double *Xbuff,double *Ybuff,double *uBuff,double *uFilt);

int gaussianFilterNonUniform(int Height,int Width,int padWidth,
                             double *Xbuff,double *Ybuff,double *uBuff,
                             double sigma2,double *uFilt);

int gaussianFilterNonUniform2(int Height,int Width,int padWidth,
                              double *Xbuff,double *Ybuff,double *uBuff,
                              double sigma2,double *uFilt);

const double gauss3x3mask[] = 
{0.0625,0.1250,0.0625,
 0.1250,0.2500,0.1250,
 0.0625,0.1250,0.0625};

const double gauss5x5mask[] = 
{0.003663003663,0.014652014652 ,0.025641025641,0.014652014652 ,0.003663003663,
 0.014652014652,0.0586080586081,0.095238095238,0.0586080586081,0.014652014652,
 0.025641025641,0.0952380952381,0.150183150183,0.0952380952381,0.025641025641,
 0.014652014652,0.0586080586081,0.095238095238,0.0586080586081,0.014652014652, 
 0.003663003663,0.014652014652 ,0.025641025641,0.014652014652 ,0.003663003663};

const double gauss7x7mask[] = 
{0.00000067,0.00002292,0.00019117,0.00038771,0.00019117,0.00002292,0.00000067,
 0.00002292,0.00078634,0.00655965,0.01330373,0.00655965,0.00078633,0.00002292,
 0.00019117,0.00655965,0.05472157,0.11098164,0.05472157,0.00655965,0.00019117,
 0.00038771,0.01330373,0.11098164,0.22508352,0.11098164,0.01330373,0.00038771,
 0.00019117,0.00655965,0.05472157,0.11098164,0.05472157,0.00655965,0.00019117,
 0.00002292,0.00078633,0.00655965,0.01330373,0.00655965,0.00078633,0.00002292,
 0.00000067,0.00002292,0.00019117,0.00038771,0.00019117,0.00002292,0.00000067};