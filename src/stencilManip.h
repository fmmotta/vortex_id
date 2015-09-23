#define _USE_MATH_DEFINES

int uFieldToGradUopenFOAM(int Height,int Width,int type,float *x0,float *dx,
                          float *X,float *Y, float *uField, float *gField);

int UToGradUuniformTorus(int Height,int Width,int type,float *x0,float *dx,
                         float *X,float *Y, float *uField, float *gField);

int UToGradnonUnifFrame(int Height,int Width,int type,float *x0,float *dx,
                        float *X,float *Y, float *uField, float *gField);