//Cubic spline with derivatives
//Written by John G Baker at gsfc.nasa.gov
#include <stdio.h>
#include <stdlib.h>
#define false 0
#define true 1

int spline_natural=false;
int spline_iold=0;
int spline_findix(double xp, const double xs[],int n);
int verbose=0;

//Solves tridiagonal for the matched second derivatives z of the 3rd order interpoling polynomials at the sample points (knots).
//The condition that 1st derivatives are continuous at the knots becomes:
//  dx[i-1] z[i-1] + 2(dx[i]+dx[i-1])z[i] + dx[i]z[i+1] = r[i]
//with:
//  dx[i] = x[i+1]-x[i]
// r[i>0] = 6(dy[i]/D[i]-dy[i-1]/D[i])
//  dy[i] = y[i+1]-y[i]
//
//Solution part 1:
//  For 0<i<n-1 defines a tridiagonal matrix  with diagonals b[i],upperdiags c[i] and lower diags a[i]
//  with numbering starting from i=1 in each case
//  All rows are identical but first and last which depend on BCs
//
//  The iterative solution algorithm by LU decomposition is:
//  (with shorthand: bb = b[i] - a[i-1]*v[i-1]) 
//  v[1] = c[1]/b[1] 
//  s[1] = r[1]/b[1] 
//  v[i] = c[i]/bb
//  s[i] = (r[i] - a[i-1]*s[i-1])/bb
//
//  In the bulk: 1<i<n-2 the abc values come directly from the eqn at top:
//  a[i]=c[i]=x[i+1]-x[i]  
//  b[i]=2*(x[i+1]-x[i-1])
//  
//  The boundary conditions for the first and last row come from some condition which
//  eliminates z[0] and z[n-1] in the i=1 and i=n-2 versions of the continuity condition at top
//  This yields the values for b[1],c[1],b[n-2],a[n-3] in the first and last rows.
//  
void spline_construct(const double xs[],const double ys[],double zs[],int n){
  //xs and ys are passed in zs are passed out
  //they all must be same length n>2
  double v[n],s[n];
  double dxm,dx,dydxm,dydx,b,am1;
  int i=1;

  //if(verbose)printf("n=%i, natural=%i\n",n,spline_natural);
  //First step initialization:
  dxm=xs[1]-xs[0];
  dx=xs[2]-xs[1];
  dydxm= (ys[1]-ys[0]) / dxm;
  dydx = (ys[2]-ys[1]) / dx;
  b= 2.0*(dx+dxm);
  if(spline_natural){//natural spline
    //these depend on the first line matrix elements
    v[1] = dx/b;
    s[1] = 6.0*( dydx - dydxm )/b; 
  } else {//not-a-knot condition (effectively fourth deriv vanishing)
    v[1] = (dx-dxm) / (b-dxm);
    s[1] =  12.0*( dydx - dydxm )/b*dx/(b-dxm);
  }
  if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g\n" ,i,dx,dydx,v[i],s[i]);
  
  
  //forward loop
  for(i=2;i<n-2;i++){
    dxm=dx;
    dx=xs[i+1]-xs[i];
    dydxm=dydx;
    dydx = (ys[i+1]-ys[i]) / dx;
    b = 2.0*(dx+dxm) - dxm*v[i-1];
    v[i] = dx / b;
    s[i] = (6.0*( dydx - dydxm ) - s[i-1]*dxm)/ b ;
    if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g %g\n" ,i,dx,dydx,v[i],s[i],( dydx - dydxm )/dx);
  }	

  //Reverse loop inititalization depends on BC:
  i=n-2;
  dxm=dx;
  dx=xs[i+1]-xs[i];
  dydxm=dydx;
  dydx = (ys[i+1]-ys[i]) / dx;
  if(spline_natural){//natural spline
    //these depend on the last line matrix elements
    //everything is identical to the bulk in this case
    am1=dxm;
    b = 2.0*(dx+dxm) - dx*v[i-1];
  } else { //not-a-knot
    am1=(1-dx/dxm)*(dxm+dx);
    b = (2.0+dx/dxm)*(dx+dxm) - am1*v[i-1];
  }
  v[i] = dx / b;
  s[i] =  (6.0*( dydx - dydxm ) - am1*s[i-1])/b;
  if(verbose) printf( "i,dx,dydx,v,s,znum %i %g %g %g %g %g\n",i,dx,dydx,v[i],s[i],( dydx - dydxm )/dx);


  //reverse loop
  zs[n-1]=0;  //this is temporary, to init the loop, not the final value
  for(i=n-2;i>0;i--){
    zs[i]= s[i] - v[i]*zs[i+1];
    if(verbose)printf( "i,z: %i %g\n",i,zs[i]);
  }

  if(spline_natural){
    zs[0]=0;
    zs[n-1]=0;
  } else {
    zs[0]=zs[1]-(xs[1]-xs[0])/(xs[2]-xs[1])*(zs[2]-zs[1]);
    zs[n-1]=zs[n-2]-(xs[n-2]-xs[n-1])/(xs[n-3]-xs[n-2])*(zs[n-3]-zs[n-2]);
  }
  if(verbose)printf( "i,z: 0 %g\n",zs[0]);
  if(1+verbose){
    for(i=1;i<n-1;i++){
      if(verbose)printf ("i,z, test: %i %g %g\n",i,zs[i],(xs[i]-xs[i-1])*zs[i-1]+2*(xs[i+1]-xs[i-1])*zs[i]+(xs[i+1]-xs[i])*zs[i+1]-6*((ys[i+1]-ys[i])/(xs[i+1]-xs[i])-(ys[i-1]-ys[i])/(xs[i-1]-xs[i])));
    }
  }
}


//The interpolating function on the [i,i+1) interval is then:
//  q(x) = y[i]*ct + y[i+1]*t - t*ct*dx^2*Ax
//with 
//  dx = x[i+1]-x[i]
//   t = (x-x[i])/dx
//  ct = 1-t
//  zL = z[i]
//  zR = z[i+1]
//  dA = (zR-zL)/6
//  Ax = (zR+zL)/4 - (ct-t)*dA/2
//and the derivatives are:
//  q'(x) =  dy[i]/dx - dx*( (ct-t)*Ax + t*ct*dA)
// q''(x) =  2( Ax - (ct-t)*dA ) 
//q'''(x) =  6dA/dx =(zR-zL)/dx        
// Arguments:
// xp    = evaluation point
// xs,ys = must be identical to xs that was provided to spline_construct 
// zs    = output from spline_construct 
// n     = array lengths 
// q,dq1,dq2,dq3  = return pointers for the value and its derivatives.
void spline_intd3(double xp, const double xs[],const double ys[],const double zs[],int n,double *q, double *dq1, double *dq2, double *dq3){
  double dx,t,ct,tct,ctmt,zL,zR,dA,Ax,yL,yR;
  int i;
  i=spline_findix(xp,xs,n);
  //printf("spline_intd3: %g < xp = %g < %g, i=%i\n",xs[0],xp,xs[n-1], i);
  if(i<0){
    printf("spline_intd3: xp is out of range.\n");
    exit(1);
  }
  dx = xs[i+1]-xs[i];
  t  = (xp-xs[i])/dx;
  if(t<0 || t>1){//do we allow extrapolation?
    printf("spline_intd3 t=%7.4f",t);
  }
  ct = 1 - t ;
  tct = t*ct;
  ctmt=ct-t;
  zL=zs[i];
  zR=zs[i+1];
  dA = (zR-zL)/6.0;
  Ax = (zR+zL)/4.0 -ctmt*dA/2.0;
  yL=ys[i];
  yR=ys[i+1];
  //results:
  *q   = yL*ct + yR*t - tct*dx*dx*Ax;
  *dq1 = (yR-yL)/dx - dx * ( ctmt*Ax + tct*dA );
  *dq2 = 2 * ( Ax - ctmt*dA );
  *dq3 = (zR-zL)/dx;
  /*if(verbose){
    printf ("results for x: %g %g %g\n",xs[i],xp,xs[i+1]);
    printf("  q: %g %g %g\n", yL,*q,yR);
    printf("dq1: %g \n", (yR-yL)/dx);
    printf("qq2: %g %g %g\n", zL,*dq2,zR);
    }*/
}

//Vector version of the same thing (Faster through a bulky (eg python) interface
void spline_intd3_vec(const double xp[],int nxp, const double xs[],const double ys[],const double zs[],int n,double q[], double dq1[], double dq2[], double dq3[]){
  int i;
  for(i=0;i<n;i++)spline_intd3(xp[i], xs, ys, zs, n, q+i, dq1+i, dq2+i, dq3+i);
}

double spline_int(double xp, const double xs[],const double ys[],const double zs[],int n){
  double dx,t,ct,tct,ctmt,zL,zR,dA,Ax,yL,yR;
  int i;
  i=spline_findix(xp,xs,n);
  //printf("spline_intd3: %g < xp = %g < %g, i=%i\n",xs[0],xp,xs[n-1], i);
  if(i<0){
    printf("spline_intd3: xp is out of range.\n");
    exit(1);
  }
  dx = xs[i+1]-xs[i];
  t  = (xp-xs[i])/dx;
  if(t<0 || t>1){//do we allow extrapolation?
    printf("spline_int t=%7.4f",t);
  }
  ct = 1 - t ;
  tct = t*ct;
  ctmt=ct-t;
  zL=zs[i];
  zR=zs[i+1];
  dA = (zR-zL)/6.0;
  Ax = (zR+zL)/4.0 -ctmt*dA/2.0;
  yL=ys[i];
  yR=ys[i+1];
  //results:
  double q   = yL*ct + yR*t - tct*dx*dx*Ax;
  return q;
}

//Bisection search for interval xi array containing xp
//for internal use.
int spline_findix(double xp, const double xs[],int n){
  int done=false;
  int iL,iR,ic;
  double xc;
  //double xL,xR;
  //set bracketing limits
  iL=0;
  iR=n-1;
  ic=spline_iold;
  /*if(verbose){
    xL=xs[iL];
    xR=xs[iR];
    xc=xs[ic];
    printf( "Entering search for %g < %g < %g with x[%i < %i < %i ]= %g\n",xL,xp,xR,iL,ic,iR,xc);
    }*/

  //first we try near the last location, hoping to get lucky
  if(xs[ic]<=xp){
    iL=ic;
    if(xs[ic+1]>=xp){
      iL=ic;
      done=true;
    } else if(ic+1<iR && xs[ic+2]>=xp){
      iL=ic+1;
      done=true;
    } else { 
      iL=ic+2; //quick search failed
    }
  } else if(xs[ic-1]<=xp) { //try one more on the left
    iL=ic-1;
    done=true;
  } else {
    iR=ic-1;//quick search failed
  }
  if(!done){
    //bisection search
    //xR=xs[iR];
    //xL=xs[iL];
    while(iR-iL>1){
      ic=(iR+iL)/2;
      xc=xs[ic];
      if(xp>=xc){
	iL=ic;
	//xL=xc;
      }
      if(xp<xc){
	iR=ic;
	//xR=xc;
      }
    }
  }
  spline_iold=iL;
  return iL;
}

        
