
int genLOseenUniformList(double Gmin,double Gmax, double rmin,double rmax,
	                     double *xmin,double *xmax, unsigned long long int seed,
	                     int numVortex,double **parVortexOut);

int genLOseenSignUniformList(double Gmin,double Gmax, double rmin,double rmax,
                         double *xmin,double *xmax, unsigned long long int seed,
                         int numVortex,double **parVortexOut);

int genLOseenBinaryList(double Gmin,double Gmax, double rmin,double rmax,
	                    double *xmin,double *xmax, unsigned long long int seed,
	                    int numVortex,double **parVortexOut);

int genLOseenLucaList(double Gmin,double Gmax, double rmin,double rmax,
                      double *xmin,double *xmax, unsigned long long int seed,
                      int numVortex,double **parVortexOut);

int genLOseenNaryList(int numG, double *Glist, int numRc, double *Rclist,
                      double *xmin,double *xmax, unsigned long long int seed,
                      int numVortex,double **parVortexOut);