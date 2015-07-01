
int genLOseenUniformList(float Gmin,float Gmax, float rmin,float rmax,
	                     float *xmin,float *xmax, unsigned long long int seed,
	                     int numVortex,float **parVortexOut);

int genLOseenSignUniformList(float Gmin,float Gmax, float rmin,float rmax,
                         float *xmin,float *xmax, unsigned long long int seed,
                         int numVortex,float **parVortexOut);

int genLOseenBinaryList(float Gmin,float Gmax, float rmin,float rmax,
	                    float *xmin,float *xmax, unsigned long long int seed,
	                    int numVortex,float **parVortexOut);

int genLOseenLucaList(float Gmin,float Gmax, float rmin,float rmax,
                      float *xmin,float *xmax, unsigned long long int seed,
                      int numVortex,float **parVortexOut);

int genLOseenNaryList(int numG, float *Glist, int numRc, float *Rclist,
                      float *xmin,float *xmax, unsigned long long int seed,
                      int numVortex,float **parVortexOut);