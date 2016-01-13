#define _USE_MATH_DEFINES

int initZero(int Height,int Width, double **gFieldOut);

int gradUtoLamb(int Height,int Width, double *gField,double **sFieldOut);

int addSingleOseen(int nVortex,double *parVortex, double *x0, double *dx, 
                   int Height,int Width, double **gFieldOut);

int addUSingleOseen(int nVortex,double *parVortex, double *x0, double *dx, 
                    int Height,int Width, double **uFieldOut);

int addConstXYShear(double *x0, double *dx,
	                int Height,int Width, 
	                double v0y0,double **gFieldOut);

int initLambOseen2D(int nVortex,double *parVortex,
                    double *x0, double *dx, int Height,int Width,
                    double **sFieldOut);

int initOseenShear2D(int nVortex,double *parVortex,
                     double *x0, double *dx, int Height,int Width,
                     double v0y0, double **sFieldOut);

int addOseen2ndGrad(int nVortex,double *parVortex, double *x0, double *dx, 
                    int Height,int Width, double **gFieldOut);

int s2ndGradUtoLamb(int nVortex,double *parVortex, double *x0, double *dx,
                    int Height,int Width, double *gField,double *sField);

int s2ndGradUtoLambShear(int nVortex,double *parVortex, double *x0, 
                         double *dx,int Height,int Width, double *gField,
                         double *sField,double v0y0);

int sndGradUwFieldToLamb(int Height,int Width,double *gField,double *wField,
                         double *sField);