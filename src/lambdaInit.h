#define _USE_MATH_DEFINES

int initZero(int Height,int Width, float **gFieldOut);

int gradUtoLamb(int Height,int Width, float *gField,float **sFieldOut);

int addSingleOseen(int nVortex,float *parVortex, float *x0, float *dx, 
                   int Height,int Width, float **gFieldOut);

int addConstXYShear(float *x0, float *dx,
	                int Height,int Width, 
	                float v0y0,float **gFieldOut);

int initLambOseen2D(int nVortex,float *parVortex,
                    float *x0, float *dx, int Height,int Width,
                    float **sFieldOut);

int initOseenShear2D(int nVortex,float *parVortex,
                     float *x0, float *dx, int Height,int Width,
                     float v0y0, float **sFieldOut);

int addOseen2ndGrad(int nVortex,float *parVortex, float *x0, float *dx, 
                    int Height,int Width, float **gFieldOut);

int s2ndGradUtoLamb(int nVortex,float *parVortex, float *x0, float *dx,
                    int Height,int Width, float *gField,float *sField);

int sndGradUwFieldToLamb(int Height,int Width,float *gField,float *wField,
                         float *sField);