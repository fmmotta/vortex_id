
#define _USE_MATH_DEFINES
#define dbg_prints(label,i) printf(label" debug %d\n",(i));
#define dbg_printp(label,i,j) printf(label " debug (%d,%d)\n",(i),(j));

#define bound_check(a,b,w,h) ((a)>=0)&&((a)<(h))&&((b)>=0)&&((b)<(w))


int check_neighbours(int i,int j,int *label,int Width,int Height,
                    int *nbList);


int findEq(int element,int *eqList,int pop);

int checkEqClass(int eqClass[][NumCls],int eqPop[],int counter);

int floodFill(float *sField,int Width,int Height,int *label);