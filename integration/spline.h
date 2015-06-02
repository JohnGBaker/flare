//Cubic spline with derivatives
//Written by John G Baker at gsfc.nasa.gov
extern int spline_verbose;
extern int spline_natural;
extern int spline_iold;

void spline_construct(const double xs[],const double ys[],double zs[],int n);
void spline_intd3(double xp, const double xs[],const double ys[],const double zs[],int n,double *q, double *dq1, double *dq2, double *dq3);
void spline_intd3_vec(const double xp[],int nxp, const double xs[],const double ys[],const double zs[],int n,double q[], double dq1[], double dq2[], double dq3[]);
double spline_int(double xp, const double xs[],const double ys[],const double zs[],int n);

        
