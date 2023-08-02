#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "plot.h"

float themin(int n, float * x){
  float m = x[0];
  for (int i=0;i<n;i++) if (m>x[i]) m=x[i];
  return(m);
}

float themax(int n, float * x){
  float m = x[0];
  for (int i=0;i<n;i++) if (m<x[i]) m=x[i];
  return(m);
}

// The min minus 5%
float themin5(int n, float * x){
  float m = themin(n,x)- 0.05*(themax(n,x)-themin(n,x));
  return(m);
}

// The max plus 5%
float themax5(int n, float * x){
  float m = themax(n,x)+0.05*(themax(n,x)-themin(n,x));
  return(m);
}

int closest(int n, float x, float y, float * xx, float * yy, float dx, float dy){
  int smallest=-1;
  float ds=0.0;
  for (int i=0;i<n;i++){
    float xdist = fabs(xx[i] - x)/dx;
    float ydist = fabs(yy[i] - y)/dy;
    float test = sqrt(xdist*xdist+ydist*ydist);
    if (smallest<0.0) {
      smallest = i;
      ds = test;
    }
    if (test<ds){
      smallest = i;
      ds = test;
    }
  }
  return(smallest);
}

float fstrideabove(float * data, size_t nchan, size_t nsubint, size_t chan, float threshold){
  int count=0;
  for (size_t i=0;i<nsubint;i++)
    if (data[i*nchan+chan]>threshold) count++;
  return(float(count)/float(nsubint));
}

float medianstride(float * data, size_t nchan, size_t nsubint, size_t chan){
  if (nsubint==1) return(data[chan]);
  float * datapts = new float[nsubint];
  for (size_t i=0;i<nsubint;i++) datapts[i]=data[chan+nchan*i];  
  return(give_median((int)nsubint,datapts));
}

float fractionstride(float * data, size_t nchan, size_t nsubint, size_t chan, float fraction){
  if (nsubint==1) return(data[chan]);
  float * datapts = new float[nsubint];
  for (size_t i=0;i<nsubint;i++) datapts[i]=data[chan+nchan*i];  
  return(give_fraction((int)nsubint,datapts,fraction));
}

float madstride(float * data, size_t nchan, size_t nsubint, size_t chan){
  if (nsubint==1) return(0.0);
  float * datapts = new float[nsubint];
  for (size_t i=0;i<nsubint;i++) datapts[i]=data[chan+nchan*i];  
  return(mad((int)nsubint,datapts));
}

/*local_medians.c, a program to calculate the running median of a 1d data array
**
**ahughes Feb 97
**
*/

float mad(int npts, float * data){
  float the_median = give_median(npts, data);
  float * temp = new float[npts];
  if (temp!=NULL){
    for (int i=0;i<npts;i++) temp[i]=fabs(data[i]-the_median);
    return(give_median(npts,temp));
  }else{
    fprintf(stderr,"Cannot allocate %zu bytes in mad routine\n",
	    sizeof(float)*npts);
  }
  return(-1);
}

int local_medians(int nsub, int nlocal, float * newdata, float * median){

/* given the input array newdata, this routine fills in the median array with
the "local" median value, which is computed from nlocal points at each point. */

  // do middle points
  for (int i=nlocal;i<=nsub;i++) 
    median[i-nlocal/2-1]=give_median(nlocal,&newdata[i-nlocal]);
    
  // Set the edge values to the median nlocal/2 away from the edge.
  for (int i=0;i<nlocal/2;i++) median[i] = median[nlocal/2];
  for (int i = nsub-nlocal/2;i<nsub; i++) median[i]=median[nsub-1-nlocal/2];
  return 0;

} 

float give_median(int npts, float * newdata){

// returns median of input array

  float * work_array = new float[npts];
  float median;
 
  // copy into working array so that original data is not overwritten
  for (int i=0;i<npts;i++) work_array[i]=newdata[i];
  sort(npts,work_array);
  median=work_array[npts/2];
  delete [] work_array;

  return(median);
}

float give_fraction(int npts, float * newdata, float f){
// returns value of point a fraction f into the sorted input array
  float * work_array = new float[npts];
  float thepoint;
  // copy into working array so that original data is not overwritten
  for (int i=0;i<npts;i++) work_array[i]=newdata[i];
  sort(npts,work_array);
  thepoint=work_array[int(npts*f)];
  delete [] work_array;
  return(thepoint);
}

//
// routine to take an array of values, nsubint * nchan and
// determine the lowest value above some fraction of the array
// in any channel and use it to define a rough "bandpass"
// then if there is significant deviation ignore it when
// computing the ideal bandpass
//

void fixbandpass(float * vals, size_t nchan, size_t nsubint, int nlocal, float tol,
	      float fraction, float * ideal_bandpass){
  //first of all determine the rough bandpass in each channel by
  //find the point below which are fraction of the points.

  //fprintf(stderr,"nchan %zu nsubint %zu nlocal %d tol %f fraction %f\n",
  //	  nchan, nsubint, nlocal, tol, fraction);

  if (fraction<0 || fraction>1){
    fprintf(stderr,"fraction (%f) in bandpass makes no sense\n",fraction);
    exit(-1);
  }

  float * mads = new float[nchan];
  for (size_t i=0;i<nchan;i++){
    ideal_bandpass[i]=fractionstride(vals, nchan, nsubint, i, fraction);
    mads[i]=madstride(vals,nchan,nsubint,i);
    //fprintf(stderr,"%zu ideal %f mads %f\n",i,ideal_bandpass[i],mads[i]);
  }

  float * locals = new float[nchan];
  float * weights = new float[nchan];
  local_medians(nchan,nlocal,ideal_bandpass,locals);
  // now look through locals, and if mad> (1.0+tol)local_median delete
  // the value and replace by interpolation
  for (size_t i=0;i<nchan;i++){
    if (mads[i]>tol*locals[i]) weights[i]=0.0; else weights[i]=1.0;
    //fprintf(stderr,"i %zu local_median %f value %f mad %f weights %f\n",i,
    //	    locals[i],ideal_bandpass[i],mads[i],weights[i]);
  }
  size_t ifirst_good=0;
  size_t ifinal_good=nchan;
  //fprintf(stderr,"Finding the first good point\n");
  for (size_t i=0;i<nchan;i++){
    if (weights[i]==1.0) {
      ifirst_good=i;
      break;
    }
  }
  //fprintf(stderr,"Found the first good point %zu\n",ifirst_good);
  //fprintf(stderr,"Finding the last good point\n");
  for (size_t i=nchan-1;i>=0 && ifinal_good>nchan-1;i--)
    if (weights[i]==1.0) ifinal_good=i;
  //fprintf(stderr,"Found the last good point %zu\n",ifinal_good);
  // Do all the points until there is a good one
  //fprintf(stderr,"Filling in the edges first_good %zu final_good %zu\n",ifirst_good,ifinal_good);
  if (ifirst_good<0 || ifinal_good>=nchan) exit(-1);
  float first_good=locals[ifirst_good];
  float last_good=locals[ifinal_good];
  for (size_t i=0;i<ifirst_good;i++) ideal_bandpass[i] = first_good;
  for (size_t i=nchan-1;i>ifinal_good;i--) ideal_bandpass[i] = last_good;
  size_t ilast_good;
  size_t inext_good;
  //fprintf(stderr,"interpolating where necessary\n");
  for (size_t i=ifirst_good;i<=ifinal_good;i++){
    if (weights[i]==1.0) {
      ilast_good=i;
      ideal_bandpass[i]=locals[i];
    }
    if (weights[i]==0.0) {
      // find the next good point.
      for (size_t j=i+1;j<=ifinal_good;j++){
	if (weights[j]==1.0){
	  inext_good=j;
	  break;
	}
      }
      ideal_bandpass[i]=locals[ilast_good]+
	(locals[inext_good]-locals[ilast_good])/float(inext_good-ilast_good)*
	float(i-ilast_good);
    }
  }
  /*
  for (size_t i=0;i<nchan;i++){
    printf("bandpass %zu %f %f %f\n",i,ideal_bandpass[i],locals[i],weights[i]);
  }
  */
}

// classic fortran-style coding found from somewhere.
// I assume it is some sort of quicksort?
void sort(int n,float *ra)
{
  int l,j,ir,i;
  float rra;

  if(n<2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
     if (l > 1)	
	rra=ra[--l-1];
     else {
	rra=ra[ir-1];
	ra[ir-1]=ra[0];
	if (--ir == 0) {
		ra[0]=rra;
		return;
		}
	}
     i=l;
     j=l<<1;
     while (j <= ir) {
	if (j < ir && ra[j-1] < ra[j]) j++;
	if (rra < ra[j-1]) {
		ra[i-1]=ra[j-1];
		i=j;
		j = i<<1;
		}
	else j=ir+1;
        }
     ra[i-1]=rra;
     }
     return;
}
