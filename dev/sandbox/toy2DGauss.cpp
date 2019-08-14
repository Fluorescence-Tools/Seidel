#include <cstdio>

typedef struct {
  int length;
  double data[1];
} LVDoubleArray;

typedef struct {
	LVDoubleArray** subimage;
	int osize;
	LVDoubleArray** M;
} MGParam;

double fit2DGaussian(double* vars, MGParam* p)
{
	// vars = [x0 y0 A sigma_x sigma_y bg info wi_nowi free_fixed_bg]
    vars[8] = 42;
    LVDoubleArray *subimage = *(p->subimage);
    LVDoubleArray *M = *(p->M);
    int length = subimage->length;
    int i;
    double *imdat = subimage->data, *mdat = M->data;
    printf("size of mdat is %i", sizeof(&mdat));

    for ( i=0; i < 400; i++){
        mdat[i] = imdat[i];
//        printf("value of imdat at position %i  is %f\n",i, mdat[i]);
    };
    printf("succesfull loop completion \n");
    
  return 1;
}