
int fprintfRunParamSigned(FILE *dadosgen,long long int seed,double x0[],
                         double xf[],double dx[],double Gmin, double Gmax,
                         double rmin, double rmax, double xmin[], double xmax[], 
                         double v0y0);

int histoIncVortex(int nVortex, double *parVortex,
                   gsl_histogram *iG, gsl_histogram *iRc,
                   gsl_histogram *ia, gsl_histogram *ib);

int fprintVortex(FILE *dadosout, int run,int nVortex, double *vCatalog);

int fprintsField(FILE *dadosout,double *x0,double *dx,
                 int Width, int Height, double *sField);

int fprintLabels(FILE *dadosout,double *x0,double *dx,
                 int Width, int Height, int *label);

int genVortices(int genType,long long int seed, double xmin[],double xmax[], 
                int nFixVortex, double **parVortex,
                double Gmin,double Gmax,double rmin,double rmax,
                double numG,double numRc, double *Glist,double *Rclist);

int calcScalarField(int runType,int Height,int Width,double x0[],double dx[],
                    int nVortex,double *parVortex,double *gField,double v0y0,
                    double *sField);

int vortexReconstruction(int runType,int Height, int Width, int nCnect, 
                         double x0[],double dx[],double *sField, 
                         double *gField,int *label,double **vCatalog);