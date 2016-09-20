void vortexAdaptiveQuickSort(double *v,int nCnect,int size,
                             int (*cmp)(const double*,const double*));

int vortexExtFromVortCurv(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut);


int vortexExtFromSwirlStr(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut);


int extractSecondMoment(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,
                        double *vCatalog,double *vortSndMomMatrix);


int extract012Momentsw2(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,
                        double *vCatalog,double *vortSndMomMatrix,
                        double *avgGradU);

int extractAvgBkgVort(int Height,int Width,double *X,double *Y,
                      int nCnect,int *label,double *wBkg,
                      double *wBkgGammaOut);

int extLambOseenParams(int Height,int Width, int nCnect,double *X,double *Y,
                       double *sField,double *gField,int *label,
                       double *vCatalog);

int extVortexVelocity(int Height,int Width, int nCnect,double *X,double *Y,
                      double *uField,double *sField,double *gField,int *label,
                      double *uVort);