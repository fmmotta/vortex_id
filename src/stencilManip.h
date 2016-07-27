#define _USE_MATH_DEFINES

int UToGradUuniformTorus(int Height,int Width,int type,double *x0,double *dx,
                         double *X,double *Y, double *uField, double *gField);

int UToGradUnonUnifFrame(int Height,int Width,int type,double *x0,double *dx,
                         double *X,double *Y, double *uField, double *gField);

int UToGradUnUnifHalfFrame(int Height,int Width,int type,double *x0,double *dx,
                           double *X,double *Y, double *uField, double *gField);

int UToGradUnUnifFullFrame(int Height,int Width,int type,double *x0,double *dx,
                           double *X,double *Y, double *uField, double *gField);