
void vortexAdaptiveQuickSort(double *v,int nCnect,int size,
                             int (*cmp)(const double*,const double*));

int vortexExtFromVortCurv(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut);


int vortexExtFromSwirlStr(int Height,int Width, int nCnect,double *X,double *Y,
                          double *sField,double *gField,int *label,
                          double **vCatalogOut);


int extractSecondMoment(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,double *vCatalog,double *vortSndMomMatrix);


int extract012Momentsw2(int Height,int Width, int nCnect,double *X,double *Y,
                        double *sField,double *gField,int *label,double *vCatalog,
                        double *vortSndMomMatrix,double *avgGradU);