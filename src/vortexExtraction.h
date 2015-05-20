
int vortexExtraction(int Height,int Width, int nCnect,
	                 float *x0, float *dx,float *sField,
	                 float *gField,int *label,float **vCatalogOut);

void vortexQuickSort(float *v,int nCnect,
                     int cmp(const float*,const float*));

int vortexExtRecursive(int Height,int Width, int nCnect,
                       float *x0, float *dx,int **eqClass,
                       float *sField,float *gField,int *label,
                       float *rCatalogOut);

int lesserCirculation(const float *v,const float *p);

int greaterCirculation(const float *v,const float *p);

int lesserVorticity(const float *v,const float *p);

int greaterVorticity(const float *v,const float *p);

int lesserRadius(const float *v,const float *p);

int greaterRadius(const float *v,const float *p);