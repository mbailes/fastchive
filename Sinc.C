/***************************************************************************
 *
 *   Copyright (C) 2005 by Russell Edwards
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <fftw3.h>
#include <stdio.h>

extern "C" {
  int ffft_(float * data, int * npts, int * isign, int * ireal);
};

// Real to complex fft of data into spectrum spec
void frc1d(int npts, float * spec, float * data){
  //  for (int i=0;i<npts;i++) spec[i]=data[i];
  //printf("Calling Joe Taylor routine with npts = %d\n",npts);
  //int one=1;
  //ffft_(spec,&npts,&one,&one);
  fftwf_plan p;
  //  printf("frc1d Creating plan for %d points\n",npts);
  p = fftwf_plan_dft_r2c_1d(npts, data, (fftwf_complex *) spec, FFTW_ESTIMATE);
  //printf("frc1d Executing FFT for %d points with FFTW3\n",npts);
  //for (int i=0;i<npts;i++) printf("frc1d %d %f %f\n",i,data[i],spec[i]);
  fftwf_execute(p); 
  //printf("frc1d Executed FFT for %d points with FFTW3\n",npts);
}

void bcr1d(int npts, float * data, float * spec){
  fftwf_plan p;
  p = fftwf_plan_dft_c2r_1d(npts, (fftwf_complex *) data, spec, FFTW_ESTIMATE);
  fftwf_execute(p);
}

#include <stdio.h>
#include <algorithm>
#include <math.h>

using namespace std;

// redwards --- code for finding the phase shift w.r.t. a template profile,
// using sinc interpolation of the cross correlation function

float sinc_interp(float *f, double x, int n)
{
  while (x < 0.0)
    x+= n;
  while (x >= n)
    x-= n;

  if (x == floor(x))
    return f[(int)floor(x)];
 
  // uses technique of Schanze, IEEE Trans Sig proc 42(6) 1502-1503 1995
  // for even n,
  //    f(x) = sin(pi.x)/pi sum_{i=-L}^{M-1} f[i] -1^i cot(pi.(x-i)/n)
  // where L+M=N.
  // We can (apparently?) choose L=0 to have M=N

  if (n%2)
  {
    fprintf(stderr, "Error, odd numbers of bins not implemented in sinc_interp!!\n"); // possible but uglier
    exit(1);
  }

  int i;
  double result = 0.0;
  for (int i=0; i < n; i++)
  {
    result += f[i] * 1.0/tan(M_PI * (x-i)/n);
    i++;
    result -= f[i] * 1.0/tan(M_PI * (x-i)/n);
  }

//  return result * sin(M_PI*x)/n;  breaks for near-integer x
  return result * sin(M_PI*(x - 2.0*floor(0.5*x)))/n;
}

static float
sinc_interp_second_derivative(float *f, double x, int n)
{
  while (x < 0.0)
    x+= n;
  while (x >= n)
    x-= n;

  if (n%2)
  {
    fprintf(stderr, "Error, odd numbers of bins not implemented in sinc_interp!!\n"); // possible but uglier
    exit(1);
  }

  int i0 = (int)floor(x+0.5);
  double delta = x - (double)i0;
  int i;

#if 1
  if (fabs(delta) < 1.0e-4)
  {

    // The equations below explode when x is nearly integer.
    // In that case, several large terms cancel and one ends up with,
    // for the peak contribution:
    double result =  (i0%2 ? -1.0: 1.0) * -f[i0] * (M_PI*M_PI/3.0 * (1.0 + 2.0/(n*n))
      + (M_PI*M_PI*(1.0+2.0*M_PI*M_PI/(3.0*n*n))) * delta * delta);
    // for the remaining terms, the following is a good appoximation:
    for (i=0; i < n; i++)
    {
      if (i!=i0)
      {
	double csc = 1.0/sin(M_PI*(x-i)/n);
	result +=  (i%2 ? -1.0: 1.0) * -f[i]*M_PI*M_PI 
	  * 2.0 * csc*csc / (1.0*n*n);
     }
    } 
    result *= (i0%2 ? -1.0: 1.0);
    return result;
//     fprintf(stderr, "Result 1 = %lg\n", result);
    //       exit(1);
  } 
#endif

  // Second Derivative of sinc_interp(), obtained using Mathematica
  double result = 0.0;
  double sinpix = sin(M_PI*(x - 2.0*floor(0.5*x)));
  double cospix = cos(M_PI*(x - 2.0*floor(0.5*x)));
  for (i=0; i < n; i++)
  {
    double csc = 1.0/sin(M_PI*(x-i)/n);
    double cscsq = csc*csc;
    double cot = 1.0/tan(M_PI*(x-i)/n);
    
    result += (i%2 ? -1.0: 1.0) * -f[i]*M_PI*M_PI * 
      (2.0*n*cospix*cscsq+cot*(1.0*n*n-2.0*cscsq)*sinpix)/(1.0*n*n*n);

//   if (fabs(delta) < 1.0e-6 && i0==41)
//     printf("%d %lg APPROX2\n", i, (i%2 ? -1.0: 1.0) * -M_PI*M_PI * 
// 	   (2.0*n*cospix*cscsq+cot*(1.0*n*n-2.0*cscsq)*sinpix)/(1.0*n*n*n));
  }
//   if (fabs(delta) < 1.0e-6 && i0==40)
//   {
//     fprintf(stderr, "Result 2 = %lg\n", result);
//     double h = 0.05;
//     fprintf(stderr, "Numerical: %lg\n",
// 	    (sinc_interp(f, x-h, n) 
// 	     - 2.0*sinc_interp(f, x, n)
// 	     +sinc_interp(f, x+h, n)) / (h*h));
//     exit(1); 

//  }
  return result ;
}

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void find_peak(float *f, unsigned n, 
	  double *xmax, float *ymax)
{
  // find peak bin
  unsigned imax = std::max_element(f, f+n)-f;

  // Rest of code is hacked version of Numerical Recipes' "brent".
  double ax = imax-1.0;
  double bx = imax;
  double cx = imax+1.0; 
  double tol = 1.0e-8; // >> sqrt of double precision
  int iter;
  double a,b,d=0,etemp,p,q,r,tol1,tol2,u,v,w,x,xm;
  float fu,fv,fw,fx;
  double e=0.0;
  const int ITMAX = 100;
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0e-10;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=-sinc_interp(f, x, n);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) 
    {
      *xmax = x;
      *ymax = -fx;
      return ;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=-sinc_interp(f, u, n);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }
  fprintf(stderr, "Internal error!!\n");
  exit(1);
}

double get_shift (int nbin, float * stdamps, float * prfamps, float * shiftval,
		  float * stddev) 
{
  int nbin_std = nbin;
  int nbin_obs = nbin;
  double mismatch_shift=0.0;
  int nby2=nbin/2;

  // compute the cross-correlation
  // Note, in case of number of bins mismatch, we compute the full FFT of
  // each and only use those coefficients they have in common
  float * obs_spec = new float[nbin+2];
  float * std_spec = new float[nbin+2];
  float * ccf_spec = new float[nbin+2];

  float zero = 0.0;
  float *ccf = new float [nbin+2];

  //  printf("Calling FFTs with nbin = %d\n",nbin_obs);
  //for (int i=0;i<nbin_obs;i++) printf("get_shift %d %f %f\n",i,prfamps[i],stdamps[i]);

  frc1d (nbin_obs, (float*)obs_spec, prfamps);
  //printf("Called FFT of profile with nbin = %d\n",nbin);
  frc1d (nbin_std, (float*)std_spec, stdamps);
  //printf("Called FFTs of std profile with nbin = %d\n",nbin);

  // complex multiplication
  // (a+bi) * (c + di)
  // real = ac -bd
  // imag = ad + bc
  // but we want (a+bi) * (c-di) = ac+bd, bc-ad

  for (int i=1; i < nby2; i++){
    // real ac + bd
    ccf_spec[2*i] = obs_spec[2*i]*std_spec[2*i]+obs_spec[2*i+1]*std_spec[2*i+1];
    // imag bc-ad
    ccf_spec[2*i+1] = obs_spec[2*i+1]*std_spec[2*i]-obs_spec[2*i]*std_spec[2*i+1];
  }
  //fprintf(stderr, "DC=%f \n", obs_spec[0]);
  ccf_spec[0] = zero; // ignore DC components
  ccf_spec[1] = zero;
  ccf_spec[nby2*2] = zero; // Ignore Nyquist, it has no phase information
  ccf_spec[nby2*2+1] = zero; // Ignore Nyquist, it has no phase information
  
  //for (int i=0;i<nbin;i++)
  //    printf("ccf_spec [%d] = %f\n",i,ccf_spec[i]);

  bcr1d (nbin, (float*)ccf_spec, ccf);

  //printf("returned from backwards fft\n");
  //for (int i=0;i<nbin;i++)
  //     printf("ccf [%d] = %f\n",i,ccf[i]);

  //fprintf(stderr, "ccf[0] = %g\n", ccf[0]);

  double maxbin;
  float peakval;
  find_peak(ccf, nbin, &maxbin, &peakval);

  peakval /= nbin*nbin; // correct for scaling in FFTs

  // Get the RMS by a means analogous to Taylor's frequency domain method.
  // We have obs = b.std + noise1, so the CCF peak, i.e. the covariance,
  //    covar = b var_std + noise2
  // and the noise variance can be got at via
  //    var_obs = b^2 var_std + var_noise1 .
  // This introduces noise in the covariance (via noise2), which will have
  // a variance of 
  //    var_noise2 = var_noise1 . var_std / N
  // Recall that 
  // so the variance of the ccf peak is
  //  (var_obs - b^2 var_std) . var std /N ... substitute b = covar / var_std
  //  = (var_obs.var_std - covar^2 ) / N

  //  variance_noise = ACF_obs(0) + b^2ACF_std(0) - 2 b CCF(tau)
  //  where b = CCF(tau) / ACF_std(0)
  double variance_obs=0.0, variance_std=0.0;
 
  for (int i=1; i < nbin/2; i++)
  {
    variance_obs += obs_spec[2*i]*obs_spec[2*i]+obs_spec[2*i+1]*obs_spec[2*i+1];
    variance_std += std_spec[2*i]*std_spec[2*i]+std_spec[2*i+1]*std_spec[2*i+1];
  }

  //   // correct for FFT scaling, plus use of only +ve frequencies
  variance_obs *= 2.0/(nbin*nbin);
  variance_std *= 2.0/(nbin*nbin);

  int nadd = nby2-1;
  double rms = (variance_obs*variance_std-peakval*peakval)/(nadd*2);

  // Arrgh for some reason this just doesn't work!! Instead copy
  // the result from chi squared minimisation, which is related to our
  // peakval by constant - 2b/var_noise1 * peakval. The delta-Chi-sq=1
  // points therefore correspond to 
  //      delta_peakval = var_noise/(2b)
  // From var_noise = var_obs - b^2 var_std, we get
  //   delta_peakval = 1/2 (var_obs/b - b var_std)
  // Substitute b = peakval / var_std:
  //  rms =-0.5*(variance_obs*variance_std/peakval - peakval);
  //  which is related to my calculation above by
  rms *= 0.5/peakval; // root 2 is fudge factor!! figure it out!!

  //  printf("rms is %f\n",rms);

  if (rms < 0.0)
  {
    fprintf(stderr, "Ooopsy\n");
  }

//   double scale = 0.5*peakval / variance_std;
//   double rms = variance_obs + scale*scale*variance_std - scale * peakval;
//   // above = var_obs + 1/4 covar^2/var_std - 1/2 covar^2/var_std
//   //       = var_obs - 1/2 covar^2/var_std


//   rms *= 1.0/nbin;

//      printf("%lg %lg %g %lg VAR\n", variance_obs, variance_std, peakval, rms);
//   printf("%lg %lg  CMP\n", 
// 	 (variance_obs*variance_std-peakval*peakval)/nadd,
// 	 rms);

 
  // now we need the second derivative of the ccf to get an idea of the 
  // error in its peak .. 
  // (note correct again for scaling in ccf)
  double second_deriv 
    = sinc_interp_second_derivative(ccf, maxbin, nbin) / (nbin*nbin);

  // the ccf in the vicinity of the peak is 
  //    ccf = ccf(maxbin) + 1/2 second_deriv (m-maxbin)^2
  // then, the change in maxbin needed to change the ccf by 1 sigma is
  //  (m-maxbin)^2 = 2. rms / second_deriv
  double sigma_maxbin =  sqrt(2.0 * rms / -second_deriv);

  //printf("rms %lg -second_deriv %lg \n", rms, -second_deriv);

  double shift = maxbin / nbin - mismatch_shift;
  double sigma_shift = sigma_maxbin / nbin;
  if (shift >=0.5)
    shift -= 1.0; 

  delete [] ccf;
  delete [] ccf_spec;
  delete [] std_spec;
  delete [] obs_spec;

  *shiftval=shift;
  *stddev=sigma_shift;

  return(0.0);  // to keep compiler happy
}

