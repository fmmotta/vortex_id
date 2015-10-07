
int vortexExtraction(int Height,int Width, int nCnect,
	                 float *x0, float *dx,float *sField,
	                 float *gField,int *label,float **vCatalogOut);

int vortexExtSimple(int Height,int Width,float *x0, float *dx,
                    int **eqClass,float *sField,float *gField,int *label,
                    float threshold,int *nCnectOut,float **vCatalogOut);

void vortexQuickSort(float *v,int nCnect,
                     int cmp(const float*,const float*));

int vortexExtRecursive(int Height,int Width,float *x0, float *dx,int **eqClass,
                       float *sField,float *gField,int *label, float threshold, 
                       float *vCatalog, int *rCnectOut,float **rCatalogOut);

int applySwirlingStrengthThreshold(int Height,int Width,float *sField,
                                   float theta);

int vortexExtTreshold(int Height,int Width, int nCnect,float theta,
                      float *x0, float *dx,float *sField,
                      float *gField,int *label,float **vCatalogOut);

int vortexExt2ndSwirl(int Height,int Width, int nCnect,
                      float *x0, float *dx,float *sField,
                      float *gField,int *label,float **vCatalogOut);

int lesserCirculation(const float *v,const float *p);

int greaterCirculation(const float *v,const float *p);

int greaterAbsCirculation(const float *v,const float *p);

int lesserVorticity(const float *v,const float *p);

int greaterVorticity(const float *v,const float *p);

int greaterAbsVorticity(const float *v,const float *p);

int lesserRadius(const float *v,const float *p);

int greaterRadius(const float *v,const float *p);