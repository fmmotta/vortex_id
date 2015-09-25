#define _USE_MATH_DEFINES

int UToGradUuniformTorus(int Height,int Width,int type,float *x0,float *dx,
                         float *X,float *Y, float *uField, float *gField);

int UToGradUnonUnifFrame(int Height,int Width,int type,float *x0,float *dx,
                         float *X,float *Y, float *uField, float *gField);

int UToGradUnUnifHalfFrame(int Height,int Width,int type,float *x0,float *dx,
                           float *X,float *Y, float *uField, float *gField);

int UToGradUnUnifFullFrame(int Height,int Width,int type,float *x0,float *dx,
                           float *X,float *Y, float *uField, float *gField);