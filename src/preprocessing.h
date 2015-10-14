
#define DEBUG_MODE false

#define id(i,j,k) (Nx*Ny*(k)+Nx*(j)+(i))

typedef struct openFoamIcoData {
    double u,v,w,p;
} openFoamIcoData;

int comp (const void * elem1, const void * elem2);

int loadAxis(FILE *nFile,int Nx,int Ny,int Nz,
             double *X,double *Y,double *Z);

int printAxis(FILE *nFile,int Nx,int Ny,int Nz,
              double  *X,double  *Y,double *Z,
              double *X2,double *Y2,double *Z2,
              const char *folder);

int loadFields(int Nx,int Ny,int Nz,FILE *uFile,FILE *pFile,
               openFoamIcoData *node);

int printYZsplitPlanes(int Nx,int Ny, int Nz,openFoamIcoData *node,
                       double *X,double *Y,double *Z,const char *folder);

int printYZcoordinates(int Nx,int Ny,int Nz,double *X,
                       double *Y,double *Z,const char *folder);