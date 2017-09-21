
#define _USE_MATH_DEFINES

#define NumCls 4096//2048

#define dbg_prints(label,i) printf(label" debug %d\n",(i));
#define dbg_printp(label,i,j) printf(label " debug (%d,%d)\n",(i),(j));

#define bound_check(a,b,w,h) ((a)>=0)&&((a)<(h))&&((b)>=0)&&((b)<(w))

int fmind(int a,int b);

int check_neighbours(int i,int j,int *label,int Width,int Height,
                    int *nbList);

int findEq(int element,int *eqList,int pop);

int checkEqClass(int eqClass[][NumCls],int eqPop[],int counter);

int floodFill(double *sField,int Width,int Height,int **eqClass,int *label);

int renameLabels(int Height,int Width,int *label);

int iterFloodFill(int Height,int Width,double *sField,int *buffer,int *label);

int iterRenameLabels(int Height,int Width,int nkeys,int *key,int *buffer,int *label);