
int vortexExtraction(int Height,int Width, int nCnect,
	                 double *x0, double *dx,double *sField,
	                 double *gField,int *label,double **vCatalogOut);

int vortexExtSimple(int Height,int Width,double *x0, double *dx,
                    int **eqClass,double *sField,double *gField,int *label,
                    double threshold,int *nCnectOut,double **vCatalogOut);

void vortexQuickSort(double *v,int nCnect,
                     int cmp(const double*,const double*));

int vortexExtRecursive(int Height,int Width,double *x0, double *dx,int **eqClass,
                       double *sField,double *gField,int *label, double threshold, 
                       double *vCatalog, int *rCnectOut,double **rCatalogOut);

int applySwirlingStrengthThreshold(int Height,int Width,double *sField,
                                   double theta);

int vortexExtTreshold(int Height,int Width, int nCnect,double theta,
                      double *x0, double *dx,double *sField,
                      double *gField,int *label,double **vCatalogOut);

int vortexExt2ndSwirl(int Height,int Width, int nCnect,
                      double *x0, double *dx,double *sField,
                      double *gField,int *label,double **vCatalogOut);

int lesserCirculation(const double *v,const double *p);

int greaterCirculation(const double *v,const double *p);

int greaterAbsCirculation(const double *v,const double *p);

int lesserVorticity(const double *v,const double *p);

int greaterVorticity(const double *v,const double *p);

int greaterAbsVorticity(const double *v,const double *p);

int lesserRadius(const double *v,const double *p);

int greaterRadius(const double *v,const double *p);