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