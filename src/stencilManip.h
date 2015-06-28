#define _USE_MATH_DEFINES

#define ux(i,j,w) (((i)*(w)+(j))*2+0)
#define uy(i,j,w) (((i)*(w)+(j))*2+1)

int checkFourStencil(int i,int j,int Width,int Height);

int updateFrom2DVelocityField(int Width,int Height,float dx,float dy,
                              float *uField,float *lambField);