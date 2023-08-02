#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include "cpgplot.h"
#include "dialog.h"
#include "plot.h"
#include "tsky.h"
#include <ctype.h>

using namespace std;

void heat();

// returns the correction in mjd for an mjd telescope pair
long double correction(long double mjd, char * telescope){
  long double dt = 0.0;
  if (strcmp(telescope,"meerkat")==0){
    if (mjd<58526.21089) return(0.0);
    if (mjd>58526.21089) dt+= (-2.4628e-5);
    if (mjd>58550.14921) dt+= (2.4630e-5);
    //this repeat is not a bug - handles two explicit corrections
    if (mjd>58550.14921) dt+= (-1.196e-6);
    if (mjd>58557.14847) dt+= (-4.785e-6);
    if (mjd>58575.95951) dt+= (5.981308411e-7);
    if (mjd>58550 &&mjd<58690) dt+= (-0.000306243);
    return(dt/86400.0); // convert to days
  }
  // otherwise no correction
  return(0.0);
}


void tolower(char * str){
  for (int i=0;i<(int)strlen(str);i++){
    if (str[i]==' ') {str[i]='\0'; break;}
    str[i]=tolower(str[i]);
  }
}

// shift between standard profile and observation prf
double get_shift(int nbin, float * std, float * prf, float * shift,
		 float * shifterror);

// Computes SEFD for MeerKAT at L-band according to Marisa
// The freq is in MHz. The SEFD includes the effect of TCMB but assumes cold sky.
float marisaSEFD(double freq, int ndish){
  return((5.71e-7 * freq*freq*freq - 1.9e-3 *freq*freq + 1.9 * freq - 113.0)/float(ndish));
}

int comma3floatscan(char * optind, float * val1, float * val2, float * val3){
  char * token = strtok(optind,",");
  int nscanned= sscanf(token,"%f",val1);
  if (nscanned==0){
    fprintf(stderr,"Error scanning float from %s\n",optind);
    return(0);
  }
  token=strtok(NULL,",");
  int nscanned2= sscanf(token,"%f",val2);
  if (nscanned2==0){
    fprintf(stderr,"Error scanning second float from %s\n",optind);
    return(1);
  }
  token=strtok(NULL,",");
  int nscanned3= sscanf(token,"%f",val3);
  if (nscanned3==0){
    fprintf(stderr,"Error scanning third float from %s\n",optind);
    return(2);
  }
  return(3);
}

int comma2floatscan(char * optind, float * val1, float * val2){
  char * token = strtok(optind,",");
  int nscanned= sscanf(token,"%f",val1);
  if (nscanned==0){
    fprintf(stderr,"Error scanning float from %s\n",optind);
    return(0);
  }
  token=strtok(NULL,",");
  int nscanned2= sscanf(token,"%f",val2);
  if (nscanned2==0){
    fprintf(stderr,"Error scanning second float from %s\n",optind);
    return(1);
  }
  return(2);
} 

int comma2intscan(char * optind, int * val1, int * val2){
  char * token = strtok(optind,",");
  int nscanned= sscanf(token,"%d",val1);
  if (nscanned==0){
    fprintf(stderr,"Error scanning int from %s\n",optind);
    return(0);
  }
  token=strtok(NULL,",");
  int nscanned2= sscanf(token,"%d",val2);
  if (nscanned2==0){
    fprintf(stderr,"Error scanning second int from %s\n",optind);
    return(1);
  }
  return(2);
} 

float modelbandpass(char * telescope, double freq, double duration){
  //  printf("bandpass telescope %s freq = %lf duration = %lf\n",telescope,freq,duration);
  // dictate the ideal bandpass as a function of frequency
  int chans[16] = {0,4,10,19,29,37,47,560,945,977,988,993,1000,1010,1020,1023};
  double bpass[16] = {1.16e-3,1.26e-3,1.63e-3,2.96e-3,5.27e-3,6.54e-3,6.81e-3,5.33e-3,
		     5.28e-3,5.48e-3,7.16e-3,7.17e-3,4.62e-3,2.78e-3,1.78e-3,
		     1.73e-3};
  double freqs[16];
  double bw=856.0;
  double frch1 = 856.0;
  
  double bpassi = bpass[0];

  for (int i=0;i<16;i++) freqs[i]=float(chans[i])*bw/1024+frch1;
  if (strncmp(telescope,"MeerKAT",7)!=0 && strncmp(telescope,"MO",2)!=0){
    fprintf(stderr,"Unrecognised telescope %s returning 1000\n",telescope);
    return(1000);
  }
  if (strncmp(telescope,"MeerKAT",7)==0){
    for (int i=0;i<15;i++) {
      if (freq>=freqs[i] && freq<=freqs[i+1]){
	bpassi = bpass[i] + (freq-freqs[i])*(bpass[i+1]-bpass[i])/(freqs[i+1]-freqs[i]);
      }
    }
    if (freq<freqs[0])bpassi = bpass[0];
    if (freq>freqs[15])bpassi = bpass[15];
  //printf("DEBUG %f %f %f\n",freq,(float)duration,bpassi/sqrt(duration/8.0));
  return(bpassi/sqrt(duration/8.0));
  }
  if (strncmp(telescope,"MO",2)==0){
    for (int i=0;i<1;i++) {
      if (freq>=800 && freq<=900){
	bpassi = 1100.0;
      }
    }
  //printf("DEBUG %f %f %f\n",freq,(float)duration,bpassi/sqrt(duration/8.0));
  return(bpassi/sqrt(duration/10.0));
  }
  return(1.0);
}

// period in seconds, ref_freq, freq in MHz, dm in normal units
// returns 0<=turns<1
float dm_shift_frac_turn(double period, double ref_freq, double freq, double dm){
  if (ref_freq==freq) return(0.0);
  float delay = 1.0/241.0 *dm*(1.0e6/freq/freq-1.0e6/ref_freq/ref_freq);
  float abs_turns = fabs(delay)/period-int(fabs(delay)/period);
  if (delay<0.0) {
    abs_turns = 1.0 - abs_turns;
  }  
  return (abs_turns);
}

class profile {
public:
  size_t nbin;
  size_t original_nbin;
  float * amps;
  double freq;
  float weight;
  profile(profile * source);
  void add(profile * p);
  void copy(profile * p);
  float residual(profile * std, int seekshift, float * rms); // removes a shifted and scale std
  void dedisperse(double cf, double dm, double folding_period);
  // but returns the S/N
  float snringate(int istart, int iwidth, float rms);
  void ascii();
  float min();
  float max();
  float rms();
  float rms(size_t istart, size_t iwidth);
  float rms(profile * std);
  int numzero();  // returns number of zero bins
  float sum();
  float sum(size_t istart, size_t n_on);
  float sum(profile * std, int val); // return the sum of the on/off pulse bins
  float mean();
  float mean(size_t istart, size_t ndat);
  float stddev(size_t istart, size_t ndat);
  float stddev();
  void scale(float factor);
  void smooth(int nbins);
  float snr_w(int nwidth, int * istart, float rms);
  float find_baseline(size_t bwidth, size_t * startbase);
  float find_baseline(profile * std);
  size_t maxbin(); // returns the maximum bin in the profile
  int maxbaseline(); // returns maximum consecutive bins in baseline
  void checkbase(); // ascii version of a profile o=on -=off
  void remove_constant(float baseline);
  void remove_baseline(size_t nbaseline);
  void remove_baseline(profile * std);
  void rotateleft(float shift, double update_freq);
  void rotateleft(float shift);
  void rotateleft(int shift);
  int shift(profile * p2);  // determine the shift in bins between profile and p2
  void subtract_amps(profile * p);  // weights ignored
  profile(size_t input_nbin, float * input_prof);
  float snr_vwauto(int * bestwidth, int * istart, float * rms, 
                 float factor, float baselinefrac, float * baseline, 
		   int minbins, int maxbins, float irms);
  float snr_vws(int * bestwidth, int * istart,
		float factor, float rms, int minbins, int maxbins);
  void bscrunch(int factor);
  void plot(int axes, int seekbounds);
  float plotresidual(profile * std, int seekshift);
};

void profile::dedisperse(double cf, double dm, double folding_period){
  float shift_turns = dm_shift_frac_turn(folding_period,cf,this->freq,dm);
  this->rotateleft(shift_turns);
  this->freq=cf;
}

int profile::numzero(){
  int counter=0;
  for (size_t i=0;i<nbin;i++) {
    if (amps[i]==0.0) counter++;
  }
  return(counter);
}

void profile::ascii(){
  for (size_t i=0;i<nbin;i++) printf("%zu %f\n",i,amps[i]);
}

float profile::mean(){
  return(this->mean((size_t)0,nbin));
}

float profile::mean(size_t istart,size_t ndat){
  float thesum = 0.0;
  for (size_t i=0;i<ndat;i++) thesum+=amps[(istart+i)%nbin];
  return(thesum/float(ndat));
}

size_t profile::maxbin(){
  size_t max_bin=0;
  float largest_bin=this->max();
  for (size_t i=0;i<nbin;i++)
    if (amps[i]==largest_bin){
      max_bin=i;
    }
  return(max_bin);
}

void profile::smooth(int nsmooth){
  float temp[nbin];
  // copy profile
  int dbin=nsmooth/2;
  for (int i=0;i<(int)nbin;i++) temp[i]=amps[i];
  
 // copy from bins +/-nsmooth/2 
  for (int i=0;i<(int)nbin;i++) {
    float sum=0;
    for (int k=0;k<nsmooth;k++){
      sum+=temp[(nbin+i+k-dbin)%nbin];
	   }
    amps[i]=sum/float(nsmooth);
  }
}

int profile::maxbaseline(){
  float lowest=this->min();
  size_t longest = 0;
  size_t stretch = 0;
  for (size_t i=0;i<2*nbin;i++){
    if (amps[i%nbin]==lowest){
      stretch++;
      //printf("-");
    } else {
      // the stretch is over, record and progress
      if (longest<stretch) longest = stretch;
      stretch=0;
      //printf("o");
    }
  }
  //  printf("\n");
  return(longest);
}

void profile::checkbase(){
  float lowest=0.0;
  size_t longest = 0;
  size_t stretch = 0;
  for (size_t i=0;i<2*nbin;i++){
    if (amps[i%nbin]==lowest){
      stretch++;
      printf("-");
    } else {
      // the stretch is over, record and progress
      if (longest<stretch) longest = stretch;
      stretch=0;
      printf("o");
    }
  }
  printf("\n");
}

// computes the flux of the profile's bins in which
//the standard is equal to val - off-pulse by default
float profile::sum(profile * std, int val=0){
  float sum=0.0;
  if (val==0){
    for (size_t i=0;i<nbin;i++){
      if (std->amps[i]==0.0) sum+=amps[i];
    }
  }
  if (val!=0){
    for (size_t i=0;i<nbin;i++){
      if (std->amps[i]!=0.0) sum+=amps[i];
    }
  }
  return(sum);
}

// computes the residual of the profile after subtracing a scaled std
// baseline is the width to use for finding the baseline
float profile::residual(profile * std, int seekshift, float * rms_baseline){
  profile * copystd = new profile(std);
  if (seekshift==0){
    // do nothing except subtract the scaled std
  } else{
    int shift = 0;
    // This takes a lot of time as of order nbin^2
    // At the end the std (copystd) is rotated to align with the profile
    if (seekshift) shift = this->shift(std);
    copystd->rotateleft(shift);
  }
  //float stdflux=copystd->sum();
  this->remove_baseline(copystd);
  *rms_baseline = this->rms(copystd);
  if (*rms_baseline==0.0) return(-1.0);
  float flux = this->sum(copystd,1);
  //this->scale(1.0/ *rms_baseline);
  //this->plot(1,1);
  //copystd->scale(flux/stdflux);  // not used now
  //this->subtract_amps(copystd);  // never used
  float snr=flux/(*rms_baseline)/sqrt(nbin-copystd->numzero());
  printf("DEBUG snr %f rms_baseline %f flux %f nbin %d numzero %d\n",
	 snr, *rms_baseline, flux, (int) nbin, copystd->numzero());
  return(snr);
}

float profile::plotresidual(profile * std, int seekshift=1){
  profile * copystd = new profile(std);
  profile * copythis = new profile(this);
  if (seekshift==0){
    // do nothing except subtract the scaled std
  } else{
    int shift = 0;
  // This takes a lot of time as of order nbin^2
    if (seekshift) shift = this->shift(std);
    copystd->rotateleft(shift);
  }
  float stdflux=copystd->sum();
  copythis->remove_baseline(std);
  //cpgsci(7);
  //copystd->plot(0,0);
  float rms_baseline = copythis->rms(std);
  copythis->scale(1.0/rms_baseline);
  cpgsci(1);
  copythis->plot(1,1);
  if (rms_baseline==0.0) return(-1.0);
  float flux = copythis->sum(copystd,1);
  copystd->scale(flux/stdflux);
  cpgsci(3);
  printf("Plotting normalised standard profile in green\n");
  copystd->plot(0,0);
  copythis->subtract_amps(copystd);
  cpgsci(2);
  // display the profile residual and compute best S/N of original
  printf("Plot residual profile in red\n");
  copythis->plot(0,0);
  float snr=flux/1.0/sqrt(nbin-copystd->numzero());
  char letters[100];
  sprintf(letters,"S/N of non-zero prf is %f\n",snr);
  cpglab(" "," ",letters);
  cpgsci(1);
  return(snr);
}

profile::profile(profile * source){
  nbin=source->nbin;
  original_nbin=source->original_nbin;
  amps=new float[nbin];
  freq=source->freq;
  weight=source->weight;
  for (size_t i=0;i<nbin;i++) amps[i]=source->amps[i];
}

void profile::remove_baseline(profile * std){
  float b=find_baseline(std);
  this->remove_constant(b);
}

void profile::remove_baseline(size_t nbaseline){
  size_t junk=0;
  float b=this->find_baseline(nbaseline,&junk);
  this->remove_constant(b);
}

// Best integer shift between profiles
int profile::shift(profile * p2){
  float sum=0.0;
  for (size_t i=0;i<nbin;i++)sum+=amps[i]*p2->amps[i];
  float best_sum=sum;
  size_t ibest_shift=0;
  for (size_t i=1;i<nbin;i++){
    sum=0.0;
    for (size_t j=0;j<nbin;j++){
      sum+=amps[j]*p2->amps[(j+i)%nbin];
    }
    if (sum>best_sum){
      best_sum=sum;
      ibest_shift=i;
    }
  }
  return(ibest_shift);
}

void profile::subtract_amps(profile * p2){
  for (size_t i=0;i<nbin;i++) amps[i]-=p2->amps[i];
}

float profile::sum(){
  return(this->sum(0,nbin));
}

float profile::sum(size_t istart,size_t n_on){
  float s=0.0;
  for (size_t i=istart;i<istart+n_on;i++) s+=amps[i%nbin];
  return(s);
}

void profile::scale(float factor){
  for (size_t i=0;i<nbin;i++) amps[i]*=factor;
}

void profile::bscrunch(int factor){
  size_t nsums = nbin/(size_t)factor;
  for (size_t i=0;i<nsums;i++){
    float sum=amps[i*(size_t)factor];
    for (size_t di=1;di<(size_t)factor;di++)
      sum+=amps[i*factor+di];
    amps[i]=sum/(float)factor;
  }
  nbin=nbin/(size_t)factor;
}

void profile::add(profile * p2){
  if (p2->weight==0.0) return;
  float sum_weight = weight+p2->weight;
  if (sum_weight==0.0) return;
  float p2weight=p2->weight;
  for (size_t i=0;i<nbin;i++){
      amps[i]=(amps[i]*weight+p2->amps[i]*p2weight)/sum_weight;
  }
  freq = (freq*weight+p2->freq*p2->weight)/sum_weight; 
  weight = sum_weight;
}

void profile::copy(profile * p2){
  for (size_t i=0;i<nbin;i++){
    amps[i]=p2->amps[i];
  }    
  freq=p2->freq;
  weight=p2->weight;
}

// assumes a pgplot window is open and viewport defined
void profile::plot(int axes, int seekbounds=1){
  float xmin = 0.0;
  float xmax = 1.0;
  float ymin = this->min();
  float ymax = this->max();
  float * xaxis = new float[nbin];
  for (int i=0;i<(int)nbin;i++) xaxis[i]=float(i)/(float(nbin-1));
  if (seekbounds) cpgswin(xmin,xmax,ymin,ymax+(ymax-ymin)*0.05);
  cpgline(nbin,xaxis,amps);
  //cpgsvp(0.15,0.85,0.15,0.85);
  if (axes) cpgbox("BCNST",0.0,0,"BC",0.0,0);
  int bestwidth=1;
  int istart=0;
  float rms=0.0, irms=0.0;
  float factor = 1.1;
  float baselinefrac=0.3;
  float baseline=0.0;
  int minbins=1;
  int maxbins=nbin*(1.0-baselinefrac);
  float snr = this->snr_vwauto(&bestwidth, &istart, &rms, 
                 factor, baselinefrac, &baseline, 
			       minbins, maxbins, irms);
  char snrinfo[100];
  float mean_flux = (this->sum()-baseline*this->nbin)/this->nbin;
  // this->min()==0 when it is a standard profile
  cpgsci(2);
  cpgmove(0.0,baseline);
  cpgdraw(1.0,baseline);
  cpgsci(1);
  if (snr<1e6 && this->min()!=0.0){
    sprintf(snrinfo,"S/N=%5.1f width %d bin strt %d",snr,bestwidth,
	    istart);
    //cpgtext(0.01,ymin+(ymax-ymin)*0.9,snrinfo);
    sprintf(snrinfo,"mean flux = %7.5f PseudoJy",mean_flux);
    //cpgtext(0.01,ymin+(ymax-ymin)*0.85,snrinfo);
    sprintf(snrinfo,"rms flux = %7.5f PseudoJy",rms);
    //cpgtext(0.01,ymin+(ymax-ymin)*0.8,snrinfo);
    cpgsls(2);
    cpgsci(7);
    cpgmove(0.0,rms+baseline);
    cpgdraw(1.0,rms+baseline);
    cpgmove(0.0,baseline-1.0*rms);
    cpgdraw(1.0,baseline-1.0*rms);
    cpgsls(1);
    cpgsci(1);
  }
}

float profile::snr_vws(int * bestwidth, int * istart,
	      float factor, float rms, int minbins, int maxbins){
  *bestwidth = minbins;   // bestwidth is initially 1.
  int testwidth = minbins;   //testwidth is initially 1
  float bestsnr = -1;   //bestsnr is initially -1
  int teststartbin = -1;   // int teststartbin records test bin
  // test_snr is set to -1;
  float test_snr = -1.0;
  int oldwidth = testwidth;
  while (testwidth<=maxbins){
    test_snr=this->snr_w(testwidth, &teststartbin, rms);
    if (test_snr>bestsnr){
      bestsnr=test_snr; *bestwidth=testwidth; *istart=teststartbin;
      //printf("snr_wvs: bestsnr %f bestwidth %d rms %f\n",test_snr,*bestwidth,rms);
    }
    testwidth=int(testwidth*factor);
    if (testwidth==oldwidth) {testwidth++;}
    oldwidth=testwidth;
  }
  return(bestsnr);
}

// Calculates the best S/N of a variable-width boxcar adjusting
// it by factor or at least one bin each time until maxbins is reached
// and assumes profile is baseline-subtracted.
// Returns the best snr, width, startbin for a given fractional baseline width
// removes the baseline prior to searching. The trial width increases
// by either 1 bin or nbin*factor-1
// pr is the profile
// nbin is the number of bins in the profile
// bestwidth is a pointer to the best width obtained for the square-pulse
// istart is the bin it starts in
// rms is the rms of the off-pulse region
// factor is the amount to increase the test number of bins by in each iteration
// baselinefrac is the fraction of the pulse to use to set the baseline
// baseline is a pointer to the amplitude of the baseline
// maxbins is the maximum width of the pulse in turns

float profile::snr_vwauto(int * bestwidth, int * istart, float * rms, 
                 float factor, float baselinefrac, float * baseline, 
		 int minbins, int maxbins, float irms){


  if (baselinefrac<=0.0 || baselinefrac>1) {
    fprintf(stderr,"profile::snr_vwauto error - baseline fraction is not 0>b>1 %f\n",baselinefrac);
    return(0.0);
  }
  // make non-destructive to input data
  profile * copyprf = new profile(this);
  size_t startbase=0;
  size_t width=int(baselinefrac*nbin);
  *baseline = copyprf->find_baseline(width,&startbase);
  copyprf->remove_constant(*baseline);
  float rms_baseline = copyprf->rms(startbase,width);
  *rms = rms_baseline;
  if (irms!=0.0) rms_baseline=irms;
  //fprintf(stderr,"profile::snr_vwauto baseline is %f rms_baseline is %f\n",*baseline,
  //	  rms_baseline);
  /*
  printf("height of baseline %f\n",baseline);
  printf("rms of baseline is %f\n",rms_baseline);
  printf("width of baseline is %d bins\n",width);
  */
  float snr=copyprf->snr_vws(bestwidth, istart, factor, rms_baseline, minbins, maxbins);
  delete copyprf;
  return(snr);
}

// Finds the best S/N of a baseline-subtracted profile of width nwidth
// and returns the start of the pulse given an input rms
float profile::snr_w(int nwidth, int * istart, float rms){
  // initial sum of nwidth bins
  if (rms==0) return(-1);
  float testsum = 0.0;
  for (int i=0;i<nwidth;i++) testsum+=amps[i];
  float bestsum=testsum;
  int ibest = 0;
  for (int i=0;i<(int)nbin;i++){
    // add a sample and subtract a sample - check for new maximum?
    testsum+=amps[(i+nwidth)%nbin];
    testsum-=amps[i];  // is mod necessary - probably not?
    if (bestsum<testsum){
      bestsum=testsum;
      ibest=i;
    }
  }
  *istart=ibest;
  return(bestsum/sqrt(nwidth)/rms);
}

// determines snr of baselined profile pr, of length nbin for width iwidth
// from starting point igate with rms rms
float profile::snringate(int iwidth, int igate, float rms){
  float sumsq=0;
  for (int i=igate;i<igate+iwidth;i++){
    int index = i%nbin;
    sumsq+=amps[index]*amps[index];
  }
  return(sqrt(sumsq/float(iwidth))/rms);
}

// shift is in turns 0<=shift<1
// rotates profile but keeps frequency the same
void profile::rotateleft(float shift){
  if (shift<0.0){
    shift = -1.0*shift;
    shift = shift - float(int(shift));
    shift = 1.0 - shift;
  }
  if (shift!=0){
    int ishift = int(shift*nbin+0.5)%(int)nbin;
    this->rotateleft(ishift);
  }
}

// rotate profile to the left by ishift bins
void profile::rotateleft(int ishift){
  if (ishift==0) return;
  int version = 2;
  if (version==1){
    if (ishift!=0){
      float temp[nbin];
      // copy profile
      for (int i=0;i<(int)nbin;i++) temp[i]=amps[i];
      for (int i=0;i<(int)nbin;i++) amps[(i+nbin-ishift)%nbin]=temp[i];
    }
  } else {
    float temp[ishift];
    // copy the soon to be obliterated values
    //    fprintf(stderr,"into temp %d\n",ishift);
    for (int i=0;i<ishift;i++) temp[i]=amps[i];
    // transfer rest
    //fprintf(stderr,"main copy %d\n",ishift);
    for (int i=0;i<(int)nbin-ishift;i++) amps[i]=amps[i+ishift];
    // tidy up at the end
    //fprintf(stderr,"tidy up %d\n",ishift);
    for (int i=0;i<ishift;i++) amps[nbin-ishift+i]=temp[i];
  }
}

// returns the mean value of the profile's bins in the baseline of std
float profile::find_baseline(profile * std){
  int nzero=std->numzero();
  float sums=this->sum(std);
  //  fprintf(stderr,"num zeroes %d sums %f\n",nzero,sums);
  if (nzero==0) {
    fprintf(stderr,"Dud standard passed to profile::findbaseline\n");
    exit(-1);
  }
  return(sums/float(nzero));
}

// returns the baseline amplitude and the starting bin of it
float profile::find_baseline(size_t bwidth, size_t * startbase){
  //determine the average of the first bwidth points
  float refsum=0.0;
  for (size_t i=0;i<bwidth;i++) refsum+=amps[i];
  float testsum = refsum;
  int testbase=0;
  // add a bin and drop a bin to find min
  for (int startbin=1;startbin<(int)nbin;startbin++){
    testsum += amps[(startbin+bwidth)%(int)nbin];
    testsum -= amps[(startbin-1)%(int)nbin];
    //printf("mean is now %f\n",testsum/float(bwidth));
    if (testsum<refsum) {
      refsum=testsum;
      testbase=startbin;
    }
  }
  *startbase = testbase;
  return(refsum/float(bwidth));
}

// finds the rms of the std's off-pulse region
float profile::rms(profile * std){
  float sumsq = 0.0;
  int icount = 0;
  for (size_t i=0;i<nbin;i++) {
    if (std->amps[i]==0.0){
      //      printf("%zu %f\n",i,this->amps[i]);
      sumsq+=pow(this->amps[i],2.0);
      icount++;
    }
  }
  if (icount>0){
    return (sqrt(sumsq/float(icount)));
  }
  else
    {
      fprintf(stderr,"profile::rms warning: no zero-bin amps in standard\n");
    }
  return(-1.0);
}

float profile::rms(size_t istart,size_t iwidth){
  float sumsq = 0.0;
  for (size_t i=0;i<iwidth;i++) sumsq+=amps[(i+istart)%nbin]*amps[(i+istart)%nbin];
  return (sqrt(sumsq/float(iwidth)));
}

float profile::rms(){
  return(this->rms(0,nbin));
}

float profile::stddev(){
  return(this->stddev(0,nbin));
}

float profile::stddev(size_t istart, size_t ndat){
  // work out the mean of the section from bin istart for ndat bins
  // careful to factor in overruns.
  float prf_mean=this->mean(istart,ndat);
  float sumsquares = 0.0;
  for (size_t i=istart;i<ndat;i++)
    sumsquares+=(amps[i%nbin]-prf_mean)*(amps[i%nbin]-prf_mean);
  return(sqrt(sumsquares/(ndat-1)));
}

void profile::remove_constant(float baseline){
  for (size_t i=0;i<nbin;i++) amps[i]-=baseline;
}

// Doesn't create the space for profiles
profile::profile(size_t input_nbin, float * input_profile){
  nbin = input_nbin;
  amps = input_profile;
  original_nbin = input_nbin;
}

float profile::min(){
  float test = amps[0];
  for (size_t i=0;i<nbin;i++){
    if (amps[i]<test) test=amps[i];
  }
  return(test);
}

float profile::max(){
  float test = amps[0];
  for (size_t i=0;i<nbin;i++){
    if (amps[i]>test) test=amps[i];
  }
  return(test);
}

// freq stored in profiles[0]->freq
class polnprofile {
public:
  size_t npol;
  size_t nbin;
  size_t original_npol;
  size_t original_nbin;
  profile ** profiles; // Max = 4
  polnprofile(size_t nbin,size_t input_npol);
  void ReWeight(double chbw, double duration);
  void set_freq(double * freq);
  void set_freq(double freq);
  void set_weight(float * w);
  void set_weight(float w);
  void copy(polnprofile * p2);
  void add(polnprofile * p2);
  void rotateleft(float shift);
  void rotateleft(float shift, double newfreq);
  void rotateleft(int ishift);
  void dedisperse(double centre_freq,double dm,double folding_period);
  void pscrunch(int nadd);
  void bscrunch(int factor);
};

void polnprofile::pscrunch(int nadd){
  if (nadd==2 && npol>1) profiles[0]->add(profiles[1]);
  profiles[0]->scale(2.0);
  npol=1;
}

void polnprofile::bscrunch(int factor){
  for (size_t i=0;i<npol;i++) profiles[i]->bscrunch(factor);
  nbin=nbin/factor;
}

void polnprofile::set_freq(double f){
  for (int i=0;i<(int)npol;i++)
    profiles[i]->freq = f;
}

void polnprofile::set_freq(double * f){
  this->set_freq(*f);
}

void polnprofile::set_weight(float *w){
  this->set_weight(*w);
}

void polnprofile::set_weight(float w){
  for (int i=0;i<(int)npol;i++)
    profiles[i]->weight = w;
}

void polnprofile::dedisperse(double cf, double dm, double folding_period){
  if (cf!=profiles[0]->freq){
    float shift_turns = dm_shift_frac_turn(folding_period,cf,profiles[0]->freq,dm);
    this->rotateleft(shift_turns,cf);
    this->set_freq(&cf);
  }
}

void polnprofile::rotateleft(float shift){
  if (shift!=0.0){
    for (size_t i=0;i<npol;i++)
      profiles[i]->rotateleft(shift);
  }
}

void polnprofile::rotateleft(int ishift){
  if (ishift!=0.0){
    for (size_t i=0;i<npol;i++)
      profiles[i]->rotateleft(ishift);
  }
}

void polnprofile::rotateleft(float shift, double update_freq){
  if (shift!=0.0){
    for (size_t i=0;i<npol;i++)
      profiles[i]->rotateleft(shift);
    this->set_freq(update_freq);
  }
}

void polnprofile::add(polnprofile * p2){
  for (size_t i=0;i<npol;i++)
    profiles[i]->add(p2->profiles[i]);
}

void polnprofile::copy(polnprofile * p2){
  for (size_t i=0;i<npol;i++)
    profiles[i]->copy(p2->profiles[i]);
}

polnprofile::polnprofile(size_t input_nbin,size_t input_npol){
  nbin = input_nbin;
  original_nbin = input_nbin;
  npol = input_npol;
  original_npol = input_npol;
  profiles = new profile * [input_npol];
}

class subint {
public:
  polnprofile ** polnprofiles;
  size_t nchan;
  size_t original_nchan;
  size_t nbin;
  size_t original_nbin;
  size_t npol;
  size_t original_npol;
  long double mjd;
  double orig_freq;   // MHz
  double orig_bw;  // MHz
  double folding_period; //spin period in seconds
  float az;
  float zen;
  float meridian();
  float meridianDeg();
  float min(size_t ipol);
  float max(size_t ipol);
  void ReWeight(double chbw, double duration);
  int nonzapped();
  float parallactic;
  double duration;
  // the functions
  double get_weighted_frequency(size_t ichan, size_t newchan);
  subint(size_t input_nchan, size_t npol);
  void copyhdr(subint * from);
  int dedisperse(double input_dm, size_t newchan);  // just rotates
  //int fadd(double input_dm);        // rotates and -F's
  int fadd(int ntoadd, double dm);  // rotates and adds ntoadd profiles together
  double get_channel_freq(size_t i);
  float get_sum_weights(size_t ichan, size_t newchan);
  float get_max_weight();
  void populate(float * data);
  void depopulate(float * data);
  void pscrunch(int nadd);
  void bscrunch(int nadd);
  void edgeZap(int nzap);
};

void subint::edgeZap(int ntrim){
  size_t current_sub=0;
  for (size_t i=ntrim;i<nchan-ntrim;i++) {
    polnprofiles[current_sub]=polnprofiles[i];
    current_sub++;
  }
  nchan = nchan - 2 * ntrim;
}

float subint::meridian(){
  // NB everything is in radians
  if (zen==0.0) return(0.0);
  float coszen = cos(zen);
  float sinzen = sin(zen);
  float sinaz = sin(az);
  float theta = atan(fabs(coszen/sinzen/sinaz));
  float absmeridian = M_PI/2-theta;
  if (az>0.0) 
    return(absmeridian);
  else
    return (0.0-absmeridian);
}

float subint::meridianDeg(){
  return (this->meridian()*180.0/M_PI);
}

//
// if the number of pols==1, set to 2*bw*t if weight>0
// or zero if weight = 0.0
//
void subint::ReWeight(double chbw, double duration){
  float factor = 1.0;
  if (this->npol==1) factor = 2.0;
  for (size_t j=0;j<nchan;j++){
    //fprintf(stderr,"Setting weight for channel %zu to %f chbw %lf duration %lf factor %f npol %d\n",
    //	    j, float(fabs(chbw*duration*factor)),chbw, duration, factor, this->npol);
    if (polnprofiles[j]->profiles[0]->weight!=0.0)
      polnprofiles[j]->set_weight(float(fabs(chbw*duration*factor)));
  }
}

void subint::copyhdr(subint * from){
  nchan=from->nchan;
  original_nchan=from->original_nchan;
  nbin=from->nbin;
  original_nbin=from->original_nbin;
  npol=from->npol;
  original_npol=from->npol;
  mjd=from->mjd;
  orig_freq=from->orig_freq;   // MHz
  orig_bw=from->orig_bw;  // MHz
  folding_period=from->folding_period; //spin period in seconds
  az=from->az;
  zen=from->zen;
  duration=from->duration;
}

int subint::nonzapped(){
  int ndata=0;
  for (size_t i=0;i<nchan;i++){
    if (polnprofiles[i]->profiles[0]->weight!=0.0) ndata++;
  }
  if (ndata==0) 
    return(0);
  else 
    return(1);
}

float subint::min(size_t ipol){
  int found_weight=0;
  float the_min;
  for (size_t i=0;i<nchan;i++){
    if (polnprofiles[i]->profiles[ipol]->weight!=0.0){
      float test=polnprofiles[i]->profiles[ipol]->min();
      if (found_weight==0) {
	the_min=test;
	found_weight=1;
      }
      if (test<the_min) the_min=test;
    }
  }
  if (found_weight) return(the_min);
  else
    return(0.0);
}

float subint::max(size_t ipol){
  int found_weight=0;
  float the_max;
  for (size_t i=0;i<nchan;i++){
    if (polnprofiles[i]->profiles[ipol]->weight!=0.0){
      float test=polnprofiles[i]->profiles[ipol]->max();
      if (found_weight==0) {
	the_max=test;
	found_weight=1;
      }
      if (test>the_max) the_max=test;
    }
  }
  if (found_weight) return(the_max);
  else
    return(0.0);
}

void subint::pscrunch(int nadd){
#pragma omp parallel for
  for (size_t i=0;i<nchan;i++)
    polnprofiles[i]->pscrunch(nadd);
  npol=1;
}

void subint::bscrunch(int factor){
#pragma omp parallel for
  for (size_t i=0;i<nchan;i++)
    polnprofiles[i]->bscrunch(factor);
  nbin=nbin/factor;
}

float subint::get_sum_weights(size_t ichan, size_t newchan){
  size_t chans2add=nchan/newchan;
  float sumweight = 0.0;
  for (size_t i=0;i<chans2add;i++){
    float w=polnprofiles[i+ichan*chans2add]->profiles[0]->weight;
    sumweight+=w;
  }
  return(sumweight);
}

float subint::get_max_weight(){
 float maxweight = polnprofiles[0]->profiles[0]->weight;
  for (size_t i=0;i<nchan;i++){
    float w=polnprofiles[i]->profiles[0]->weight;
    if (w>maxweight) maxweight=w;
  }
  return(maxweight);
}

// return the weighted frequency of the ichan'th channel 
// of newchans long
double subint::get_weighted_frequency(size_t ichan, size_t newchan){
  size_t test = nchan%(int)newchan;
  size_t chans2add=nchan/newchan;
  if (test!=0){
    fprintf(stderr,"subint::get_weighted_frequency error\n");
    fprintf(stderr,"newchan=%d must divide nchan=%d",(int)newchan,(int)nchan);
    return(0);
  }
  double sumfreqweight = 0.0;
  for (size_t i=0;i<chans2add;i++){
    float w=polnprofiles[i+ichan*chans2add]->profiles[0]->weight;
    double f=polnprofiles[i+ichan*chans2add]->profiles[0]->freq;
    //    fprintf(stderr,"ichan %d w %f f %lf\n",(int)i,w,f);
    sumfreqweight+=(double)w*f;
  }
  float sumw=this->get_sum_weights(ichan,newchan);
  //fprintf(stderr,"subint::get_weighted_frequency weight_sum %f sum_weight_freq %lf\n",sumw,sumfreqweight);
  if (sumw!=0.0) return(sumfreqweight/sumw);
  else
    return(0.0);
}
// given the mid
double subint::get_channel_freq(size_t ichan){
  return(polnprofiles[ichan]->profiles[0]->freq);
}

// dedisperse to a central frequency.
int subint::dedisperse(double input_dm, size_t newchan){
  // get the reference frequency
  //fprintf(stderr,"Weighted Centre Freq %lf MHz ",cfreq);
  // rotate each polnprofile
  size_t ntoadd = nchan/newchan;
  double cfreq;
  double f;
  float shift_turns;
#pragma omp parallel for private(cfreq,f,shift_turns)
  for (size_t i=0;i<newchan;i++){
    cfreq = this->get_weighted_frequency(i,newchan);
#pragma omp parallel for private (f,shift_turns)
    for (size_t j=0;j<ntoadd;j++){
      f = polnprofiles[i*ntoadd+j]->profiles[0]->freq;
      shift_turns = dm_shift_frac_turn(folding_period,cfreq,f,input_dm);
      polnprofiles[i*ntoadd+j]->rotateleft(shift_turns,f);
    }
  }
  return(0);
}

/* Dedisperses a subint on nbin*nchan to nearest bin
// but don't scrunch
// tested on old version of fastpdmp
void dedisperse_subint(int nbin, int nchan, float * subint, float frch1, float df,
		      float dm, float period, float * dedispersed){
    // for each channel, determine the shift
    float ref=get_ref_freq(frch1, df, nchan);
    for (int i=0;i<nchan;i++){
      float chan_freq = frch1 + df * (float)i;
      float shift=dmshift_turns(period,ref,chan_freq,dm);
      if (shift<0.0) shift = 1.0 + shift;
      int binshift = (int(shift*nbin+0.5))%nbin;
      //printf("Shift channel %d of %d is %d\n",i,nchan,binshift);
      rotate(nbin,&subint[i*nbin],binshift,&dedispersed[i*nbin]);
    }
}
*/ 

int subint::fadd(int ntoadd, double dm){
  int nadd;
  if (ntoadd==0) {
    nadd=nchan;
  } else{
    nadd=ntoadd;
  }
  if ((nchan % nadd)!=0){
    fprintf(stderr,"ntoadd %d must divide existing number of channels %zu\n",
	    nadd, nchan);
    return(-1);
  }
  size_t new_nchan = nchan/nadd;
  for (int i=0;i<(int)new_nchan;i++){
    // determine new centre frequency of each sub-band
    //fprintf(stderr,"weighted frequency for channel %d is now %lf\n",
    //	    (int)i,fcentre);
    // dedisperse the profiles in each channel
    double fcentre = this->get_weighted_frequency(i,new_nchan);
    #pragma omp parallel for
    for (int j=0;j<nadd;j++)
      polnprofiles[j+nadd*i]->dedisperse(fcentre,dm,folding_period);
    // now copy the first profile into destination
    polnprofiles[i]->copy(polnprofiles[i*nadd]);
    // now add the rest.
    for (int j=1;j<nadd;j++)
      polnprofiles[i]->add(polnprofiles[i*nadd+j]);
    // update the centre frequency of the subint and profile
    polnprofiles[i]->set_freq(&fcentre);
  }
  nchan=new_nchan;
  return(new_nchan);
}

// dangerous in the all polnprofiles are not pointing anywhere
// must be followed up by creation of the profiles.
subint::subint(size_t input_nchan, size_t input_npol){
  nchan = input_nchan;
  original_nchan = nchan;
  npol = input_npol;
  original_npol = input_npol;
  polnprofiles = new polnprofile * [input_nchan];
}

void subint::populate(float * f){
  folding_period= * (double *) &f[0];
  mjd = * (long double *) &f[2];
  az = f[6];
  zen = f[7];
  parallactic = f[8];
  duration = * (double *)&f[9];
  //  printf("period %lf duration %lf\n",folding_period,duration);
}

// Moves the subint header into contiguous RAM
void subint::depopulate(float * f){
  memcpy(&f[0],&folding_period,sizeof(double));
  //folding_period= * (double *) &f[0];
  memcpy(&f[2],&mjd,sizeof(long double));
  //mjd = * (long double *) &f[2];
  memcpy(&f[6],&az,sizeof(float));
  memcpy(&f[7],&zen,sizeof(float));
  memcpy(&f[8],&parallactic,sizeof(float));
  //duration = * (double *)&f[9];
  memcpy(&f[9],&duration,sizeof(double));
}

class archive {
public:
  size_t nbin;
  size_t npol;
  size_t nchan;
  size_t nsubint;
  size_t original_nbin;
  size_t original_npol;
  size_t original_nchan;
  size_t original_nsubint;
  float freq;
  float bw;
  double dm;
  char source[16];
  char telescope[16];
  double ra;
  double dec;
  char * filename;
  // double period;    // period of the pulsar
  // double dt;        // subint gap in seconds
  float * data;
  subint ** subints;
  archive(char * filename);                 // opens a file and reads in an archive
  archive(char * data, size_t length);     // gets the address of an archive
  archive(archive * from);                // copies an existing archive
  void refill(archive * src);             // refills an archive from src with no mallocs 
  void dedrift(double dTurnsdDay, long double tref); // rotates each subint by dTurnsdDay * deltaT
  void edgeZap(float edgefrac);
  void maxMeridianAngle(float maxang);
  void incoadd(archive * target);  // incoherently adds together across antennas
  void set_DM(double input_dm);
  void ReWeight();
  void trim(float tmin);
  //  void newPeriod(double newperiod); // rephases archive based on a new centre folding period
  float min(size_t ipol);
  float max(size_t ipol);
  double get_weighted_frequency(size_t chan, size_t isubstart, size_t nsub);
  double get_weighted_frequency(size_t chan);
  double get_weighted_frequency();
  long double get_median_mjd(size_t istart,size_t nadd);
  double get_total_duration(size_t istart,size_t nadd);
  double get_folding_period_median(size_t istart,size_t nadd);
  int showHeader(int verbose);
  int showIntHeader();
  int showWeightHeader();
  int unload();
  int unload(char * newname);
  int unload(char * root, char * extension);
  void rotateleft(float shift);
  void rotateleft(int nshift);
  void remove_baseline(size_t ipol, int nbaseline);
  void tscrunch(size_t nadd);
  void fadd(int ntoadd);
  void bscrunch(int scrunchfactor);  // must divide nbin
  void pscrunch(int plan);  // plan = 1 I, plan=2 LL+RR, default=1  
  double dedisperse(); // dedisperse to a central frequency - no scrunching
  double dedisperse(double input_dm); // dedisperse to a central frequency - no scrunching
  double dedisperse(double input_dm, double refFreq); // dedisperse to given dm and frequency
  void multiplot();
  void update_source(char * newpsrname);
  void Yplot(size_t ipol, int click);
  void Gplot(size_t ipol, int click);
  void Splot(size_t isub, size_t ichan, int click);
  void Splot(int click);
  void explore(char * header, float ** metadata, profile * std);
  void mkstd();
  void pdmp(size_t nbin, size_t nsub, size_t nchan, float dPturns, float dDMturns,
	    float wmin, float wmax, float wfactor, float basewidth, float graininess,
	    int *bestw, float * bestSNR,
	    float * bestP, float * bestDM, profile * bestprf, int * beststart);
};


// multi-threaded pdmp - scrunches archive to nbin*nsub*nchan
// then searches for appropriate slopes for the maximum S/n for
// a given width search plan given a baseline fractional width of basewidth
// graininess of 1.0 means advance by a full width of pulse between trials
// graininess = 0.5 is more conservative

float blah(){
  return(1);
}

void archive::update_source(char * newpsrname){
  if (strlen(newpsrname)<16)
    strncpy(this->source,newpsrname,16);
  else
    fprintf(stderr,"archive::update_source error source >>%s<< is too long %d chars\n",
	    newpsrname,(int)strlen(newpsrname));
}

// Performs a pdmp-like search on an archive returning the best width, SNR, Period, DM
// and a copy of the best profile and the start of the window that contains the flux.
// The search firstly scrunches in bins, subints and chans to speed things up.
// dDPturns is the number of phase turns over which the drift can occur
// dDMturns is the number of phase turns over which the DM drift can occur
// wstart,wend and wfactor help only search over a sensible width
// basew is the baseline to use for calculating the rms of the baseline
// graininess is some relation to how coarse the search is
void archive::pdmp(size_t nbint, size_t nsubt, size_t nchant, float dPturns, float dDMturns,
		   float wstart, float wend, float wfactor, float basew, float graininess,
		   int *bestw, float * bestSNR, 
		   float * bestP, float * bestDM, profile * bestprf, int * beststart){

  // First of all check we can scrunch to the correct dimensions
  // nbin first
  int verbose=1;

  if (nbint>this->nbin || (this->nbin%nbint !=0)){
    fprintf(stderr,"nbin of %zu is not appropriate for archive of nbin %zu\n",
	    nbint,this->nbin);
    return;
  }

  if (nsubt > this->nsubint){
    fprintf(stderr,"nsubint_target of %zu is not appropriate for archive of nsubint %zu\n",
	    nsubt,this->nsubint);
    return;
  }

  if (nchant>this->nchan || (this->nchan%nchant !=0)){
    fprintf(stderr,"nchan_target of %zu is not appropriate for archive of nchan %zu\n",
	    nchant,this->nchan);
    return;
  }

  //make a copy of the archive so op is non-destructive and scrunch it down

  archive * a = new archive(this);
  a->bscrunch(nbint/a->nbin);
  a->fadd(nchant/a->nchan);
  a->tscrunch(nsubt/a->nsubint);

  // Now loop over DMs and periods trying to optimise fit

  // DMsteps depends upon minimum width and range to search in turns  
  // time steps depends on time range and minimum width to search in turns
  // these are the increments in turns for each dimension

  float dshiftPturns = wstart*graininess;
  float dshiftDMturns = wstart*graininess;

  // Find the total duration and set a mid time about which to pivot
  // and assume a->freq is the centre frequency

  long double midmjd = (a->subints[a->nsubint-1]->mjd+a->subints[0]->mjd)/2.0;
  long double halfdmjd = (a->subints[a->nsubint-1]->mjd-a->subints[0]->mjd)/2.0;
  long double dPhaseBydDay=0.0;
  if (a->nsubint>1) dPhaseBydDay = dPturns/halfdmjd;   // maximum phase shift per day

  // determine the DM step - which is the dm step at the lowest frequency
  // delay = (1/nu^2-1/cf^2) * dDM * K
  // so dDm = delay/K/(1/nu^2-1/cf^2)=dshiftPturns*pfold/K/(1/nu^2-1/cf^2);

  double nulow = a->freq-fabs(a->freq-a->subints[0]->polnprofiles[0]->profiles[0]->freq);
  if (a->bw<0.0){
    nulow = a->freq-fabs(a->freq-a->subints[0]->polnprofiles[a->nchan-1]->profiles[0]->freq);
  }
  double cf = a->freq;
  float dDM = 0.0;
  if (cf>nulow && a->nchan>1){
    //      double quotient = 
      dDM = dshiftDMturns*a->subints[a->nsubint/2]->folding_period*241.0/
	(1.0e6/(nulow*nulow)-1.0e6/(cf*cf));
  }
  //  printf("dDM is %f nulow is %lf cf is %lf\n",dDM,nulow,cf);
  // create space for enough tscrunch trials to each have their own profile to aid parallelisation
  int ndmsteps = dDMturns/dshiftDMturns;
  if (a->nchan==1) ndmsteps=0;
  int nperiodsteps = dPturns/dshiftPturns;
  if (a->nsubint==1) nperiodsteps=0;
  //if (verbose) fprintf(stderr,"dm steps %d period steps %d\n",ndmsteps,nperiodsteps);
  float snrs[(2*ndmsteps+1)][(2*nperiodsteps+1)];
  float widths[(2*ndmsteps+1)][(2*nperiodsteps+1)];
  float allsnrs[(2*ndmsteps+1)*(2*nperiodsteps+1)];
  int allsnrscount=0;
  fprintf(stderr,"Ndm steps %d Nperiod steps %d\n",ndmsteps*2+1,nperiodsteps*2+1);

  // search the zero, zero archive and make spaces
  archive * dmsearch = new archive(a);
  dmsearch->dedisperse();
  dmsearch->fadd(0);
  archive * psearch = new archive(dmsearch);
  psearch->tscrunch(0);
  int bestwidth;
  int istart;
  float rms,irms=0.0;
  float baseline;
  float	bestdm;
  float	bestdrift;
  float bestsnr=psearch->subints[0]->polnprofiles[0]->profiles[0]->
    snr_vwauto(&bestwidth, &istart, &rms,
                 wfactor, basew, &baseline, 
	       int(wstart*psearch->nbin), int(wend*psearch->nbin), irms);
  bestprf = new profile(psearch->subints[0]->polnprofiles[0]->profiles[0]);
  int dd = 0;
  //  if (verbose) fprintf(stderr,"Number of Dm steps: %d number of period steps %d\n",ndmsteps,nperiodsteps);
  //  if (verbose) fprintf(stderr,"zero-zero SNR is %f\n",bestsnr);
  bestsnr=0.0;

  int bestp = 0;
  int bestd = 0;
  for (int d = -ndmsteps;d<=ndmsteps;d++){
    dmsearch->refill(a);
    double trialdm = a->dm+double(d)*dDM;
    printf("trial dm for index %d is %f\n",d,trialdm);
    dmsearch->dedisperse(trialdm);
    dmsearch->fadd(0);
    int pp = 0;
    for (int p = -nperiodsteps;p<=nperiodsteps;p++)
      {
	double pdrift = double(p)*dPhaseBydDay/double(nperiodsteps);
	psearch->refill(dmsearch);
	psearch->dedrift(pdrift,midmjd);
	psearch->tscrunch(0);
	int testwidth;
	float testsnr=psearch->subints[0]->polnprofiles[0]->profiles[0]->
	  snr_vwauto(&testwidth, &istart, &rms,
                 wfactor, basew, &baseline, 
	       int(wstart*psearch->nbin), int(wend*psearch->nbin), irms);
	// first index is dm, second is period
	snrs[d+ndmsteps][p+nperiodsteps]=testsnr;
	allsnrs[allsnrscount]=testsnr;
	allsnrscount++;
	//widths[pp+dd*(nperiodsteps*2+1)]=widthtest;
	// remember snrs for plotting and store the best profile
	if (testsnr>bestsnr){
	  *beststart=istart;
	  bestsnr=testsnr;
	  bestdm=trialdm;
	  bestdrift=pdrift;
	  bestwidth=testwidth;
	  bestprf->copy(psearch->subints[0]->polnprofiles[0]->profiles[0]);
	  bestp=p;
	  bestd=d;
	}
	pp++;
        if (ndmsteps==0)printf(" %5.2f\r",100.0*float(pp)/(nperiodsteps*2+1));
	if (ndmsteps==0)fflush(stdout);
      }
    dd++;
    if (ndmsteps!=0){
      printf(" %5.2f\r",100.0*float(dd)/(ndmsteps*2+1));
      fflush(stdout);
    }
  }
  printf("\n%s best dm %f best drift %f turns/day snr %f width %d bins\n",
	  this->filename,float(bestdm),
	  float(bestdrift),bestsnr,bestwidth);
  double delta=1.0+bestdrift*subints[a->nsubint/2]->folding_period*86400.0;
  *bestP=delta*subints[a->nsubint/2]->folding_period;
  *bestw=bestwidth;
  *bestDM=bestdm;
  *bestSNR=bestsnr;
  cpgeras();
  // plot best phase vs freq
  dmsearch->refill(a);
  dmsearch->dedisperse(bestdm);
  dmsearch->tscrunch(0);
  cpgsvp(0.1,0.3,0.65,0.95);
  dmsearch->Gplot(0,0);
  // now dedisperse it
  dmsearch->refill(a);
  dmsearch->dedisperse(bestdm);
  dmsearch->fadd(0);
  psearch->refill(dmsearch);
  psearch->dedrift(bestdrift,midmjd);
  cpgsvp(0.4,0.6,0.65,0.95);
  psearch->Yplot(0,0);
  // plot psearch to see phase vs time
  psearch->tscrunch(0);
  // plot the profile
  cpgsvp(0.1,0.9,0.1,0.25);
  psearch->Splot(0);
  // show numbers on the plot

  // plot the bullseye
  cpgsvp(0.7,0.9,0.65,0.95);
  cpgswin(0.0,float(2*nperiodsteps+1),0.0,float(2*ndmsteps+1));
  heat(); // sets the colour map for cpgimag
  float tr[6];
  tr[0]=0.0-0.5;
  tr[1]=1.0; tr[2]=0.0; tr[3]=0.0-0.5; tr[4]=0.0; tr[5]=1.0;
  cpgimag(allsnrs,nperiodsteps*2+1,ndmsteps*2+1,1,nperiodsteps*2+1,
	  1,ndmsteps*2+1,
	  themin((2*ndmsteps+1)*(2*nperiodsteps+1),allsnrs),
	  themax((2*ndmsteps+1)*(2*nperiodsteps+1),allsnrs),
	      tr);
  // Now reset 
  cpgswin(-1.0*dPturns,1.0*dPturns,-1.0*ndmsteps*dDM,1.0*ndmsteps*dDM);
  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
  cpglab("Drift","Delta DM","Bullseye");
  
  // plot SNR vs DM
  // index of array with best DM is bestd
  // An array to hold x and snr
  float * xaxis = new float[ndmsteps*2+1];
  float * plotsnr = new float[ndmsteps*2+1];
  for (int i=0;i<ndmsteps*2+1;i++){
    xaxis[i]=(i-ndmsteps)*dDM;
    plotsnr[i]=snrs[i][bestp+nperiodsteps];
  }
  cpgsvp(0.1,0.3,0.35,0.55);
  cpgswin(xaxis[0],xaxis[2*ndmsteps],0.0,*bestSNR*1.05);
  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
  cpgline(2*ndmsteps+1,xaxis,plotsnr);
  cpglab("Delta DM","SNR"," ");
  delete [] xaxis;
  delete [] plotsnr;
  xaxis = new float[nperiodsteps*2+1];
  plotsnr = new float[nperiodsteps*2+1];
  for (int i=0;i<nperiodsteps*2+1;i++){
    xaxis[i]=dPturns*(i-nperiodsteps)/nperiodsteps;
    plotsnr[i]=snrs[bestd+ndmsteps][i];
  }
  cpgsvp(0.4,0.6,0.35,0.55);
  cpgswin(xaxis[0],xaxis[2*nperiodsteps],0.0,*bestSNR*1.05);
  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
  cpgline(2*nperiodsteps+1,xaxis,plotsnr);
  cpglab("Delta P","SNR"," ");
  delete [] xaxis;
  delete [] plotsnr;
}

// only keeps subints with an absolute value of meridian angle < maxang
void archive::maxMeridianAngle(float maxang){
  size_t current_subint = 0;
  for (size_t i=0;i<nsubint;i++){
    if (fabs(subints[i]->meridianDeg())<=maxang){
      subints[current_subint]=subints[i];
      current_subint++;
    }
  }
  nsubint=current_subint;
}

// only keeps inner frequency channels 
void archive::edgeZap(float frac){
  if (frac!=0.0){
    int ntrim=0;
    if (frac>=1){
      ntrim = int(frac);
    } else {
     ntrim = int( frac* nchan);
   }
    for (size_t i=0;i<nsubint;i++){
      subints[i]->edgeZap(ntrim);
    }
    bw = bw * float(nchan - 2 * ntrim)/float(nchan);
    nchan = nchan - 2 * ntrim;
  }
}

void archive::ReWeight(){
  for (size_t i=0;i<nsubint;i++){
    //fprintf(stderr,"reweighting subint %zu\n",i);
    subints[i]->ReWeight(double(bw/nchan),subints[i]->duration);
  }
}

void archive::mkstd(){
  //pscrunch(2);
  char * psrname=strtok(source," ");
  if (npol>1) pscrunch(2);
  fadd(0);
  tscrunch(0);
  // this creates a "feature" for large duty cycle profiles
  this->subints[0]->polnprofiles[0]->profiles[0]->remove_baseline(nbin/8);
  dialog * d = new dialog();
  d->addplotregion(0.2,0.95,0.3,0.85);
  int QUIT=d->addbutton(0.02,0.02,"QUIT");
  int PLOT=d->addbutton(0.02,0.06,"PLOT");
  int BACKTOORIGIN=d->addbutton(0.5,0.06,"BACK TO ORIGIN");
  int EDGEIT=d->addbutton(0.02,0.1,"MOVE PEAK TO EDGE");
  int SMOOTH3=d->addbutton(0.02,0.14,"SMOOTH3");
  int SMOOTH7=d->addbutton(0.02,0.18,"SMOOTH7");
  int SMOOTH15=d->addbutton(0.02,0.22,"SMOOTH15");
  int SMOOTH31=d->addbutton(0.02,0.26,"SMOOTH31");
  char * buttonname = new char[strlen(psrname)+5+9];
  sprintf(buttonname,"SAVE %s.raw.std",psrname);
  int SAVESTD=d->addbutton(0.7,0.06,buttonname);
  float min=this->subints[0]->polnprofiles[0]->profiles[0]->min();
  float max=this->subints[0]->polnprofiles[0]->profiles[0]->max();
  cpgeras();
  cpgsci(1);
  d->plotregions[0].set(0.0,1.0,min,max);
  d->draw();
  cpgsci(1);
  d->plotregions[0].reset();
  this->subints[0]->polnprofiles[0]->profiles[0]->plot(1,0);
  cpglab("Pulsar Phase","Amplitude","Click on left gate and right gate of baseline");
  int button=-1;
  int plotno;
  float x,y;
  char ans;
  size_t ithbin=0;
  while ((button=d->manage(&x,&y,&ans,&plotno))!=QUIT){
    if (ans=='q') button=QUIT;
    if (button==SAVESTD){
      char * ext=new char[9];
      sprintf(ext,".raw.std");
      // no bin shall be less than zero.
      for (size_t ibin=0;ibin<nbin;ibin++) {
	if (subints[0]->polnprofiles[0]->profiles[0]->amps[ibin]<0.0) 
	  subints[0]->polnprofiles[0]->profiles[0]->amps[ibin]=0.0001*max;
      }
      int i=0;
      printf("Unloading new standard profile %s with extension %s %d\n",psrname,ext,i);
      unload(psrname,ext);
      i++;
      printf("Done! Now quitting\n");
      button=QUIT;
      break;
    }
    if (button==SMOOTH3){
      subints[0]->polnprofiles[0]->profiles[0]->smooth(3);
      button=PLOT;
    }
    if (button==SMOOTH7){
      subints[0]->polnprofiles[0]->profiles[0]->smooth(7);
      button=PLOT;
    }
    if (button==SMOOTH15){
      subints[0]->polnprofiles[0]->profiles[0]->smooth(15);
      button=PLOT;
    }
    if (button==SMOOTH31){
      subints[0]->polnprofiles[0]->profiles[0]->smooth(31);
      button=PLOT;
    }
    if (button==EDGEIT){
      ithbin=subints[0]->polnprofiles[0]->profiles[0]->maxbin();
      subints[0]->polnprofiles[0]->profiles[0]->rotateleft((int)ithbin);
      button=PLOT;
    }
    if (button==BACKTOORIGIN){
      subints[0]->polnprofiles[0]->profiles[0]->rotateleft(int(nbin-(int)ithbin));
      button=PLOT;
    }
    if (button==PLOT) {
      cpgsci(1);
      d->plotregions[0].erase();
      d->plotregions[0].reset();
      this->subints[0]->polnprofiles[0]->profiles[0]->plot(1,0);
    }
    if (plotno==0){
      int XMODE=4;
      int POSN=0;
      float XREF=x,YREF=y;
      d->plotregions[0].reset();
      cpgsci(2);
      cpgband(XMODE,POSN,XREF,YREF,&x,&y,&ans);
      // delete amps from XREF to x then replot
      for (int i=XREF*nbin;i<x*nbin;i++){
	subints[0]->polnprofiles[0]->profiles[0]->amps[i]=0.0;
	float xdat=i/(float)nbin;
	cpgslw(5);
	cpgsch(3.0);
	cpgpt1(xdat,0.0,2);
	cpgslw(1);
	cpgsch(1.0);
	usleep(1000);
      }
      cpgsci(2);
      cpgslw(5);
      cpgmove(XREF,0.0);
      cpgdraw(x,0.0);
      cpgsci(1);
      cpgslw(1);
      button=PLOT;
    }
    if (button==PLOT) {
      cpgsci(1);
      d->plotregions[0].erase();
      d->plotregions[0].reset();
      this->subints[0]->polnprofiles[0]->profiles[0]->plot(1,0);
    }
  }
}

long double archive::get_median_mjd(size_t start_sub,size_t nadd){
  // 1. If there is an odd number of subints set to middle one
  // 2. else give the mean of the two nearest modulo mean of folding period
  long double mjd_middle;
  size_t end_sub = start_sub+nadd-1;
  if (end_sub>=nsubint) end_sub=nsubint-1;
  // OK - even or odd?
  int even = (end_sub-start_sub) % 2;
  if (!even || start_sub==end_sub) return(subints[(start_sub+end_sub)/2]->mjd);
  // otherwise set to the mean of innermost two subints
  // modulo mean folding period
  if (even){
    size_t lower = (start_sub+end_sub)/2;
    size_t upper = lower+1;
    long double mjd1=subints[lower]->mjd;
    long double mjd2=subints[upper]->mjd;
    mjd_middle=(mjd1+mjd2)/2.0;
    // now compute the "whole" number of turns between these to set
    // the time to.
    double fp1=subints[lower]->folding_period;
    double fp2=subints[upper]->folding_period;
    double fpmean=(fp1+fp2)/2.0;
    size_t turns = (mjd_middle-mjd1)*86400/(long double)fpmean;
    if (turns<=1) {
      // just return mjd1
      return(mjd1);
    }
    if (turns>1){
      //compute new half number of turns
      size_t nt=turns/2;
      // the new MJD is just the turn before the middle one
      return(mjd1+((long double)nt*fpmean)/86400.0);
    }
  }
  fprintf(stderr,"Reached fatal point in archive::get_median_mjd\n");
  exit(-1);
}
double archive::get_total_duration(size_t startsub,size_t nadd){
  double total=0.0;
  for (size_t i=startsub;i<nadd && i<nsubint;i++)total+=subints[i]->duration;
  return(total);
}
// trim the last archive if length<mintime
void archive::trim(float mintime){
  int newnsub=0;
  if (subints[nsubint-1]->duration<mintime)nsubint=nsubint-1;
}

double archive::get_folding_period_median(size_t start_sub,size_t nadd){
  // 1. If there is an odd number of subints set to middle one
  // 2. else give the mean of the two nearest
  double fpmean;
  size_t end_sub = start_sub+nadd-1;
  if (end_sub>=nsubint) end_sub=nsubint-1;
  // OK - even or odd?
  int even = (end_sub-start_sub) % 2;
  if (!even || start_sub==end_sub) fpmean = subints[(start_sub+end_sub)/2]->folding_period;
  // otherwise set to the mean of innermost two subints
  if (even){
    size_t lower = (start_sub+end_sub)/2;
    size_t upper = lower+1;
    double fp1= subints[lower]->folding_period;
    double fp2=subints[upper]->folding_period;
      fpmean=(fp1+fp2)/2.0;
  }
  return(fpmean);
}

int parse(char * header,const char * delimiter){
  char * token;
  int count = 0;
  char * copyheader = new char[strlen(header)+1];
  strcpy(copyheader,header);
  token = strtok(copyheader,delimiter);
  while (token!=NULL) {
    count++;
    printf("parsed token %d %s\n",count,token);
    token = strtok(NULL, delimiter);
  }
  return(count);
}

// uses dialog to explore arrays with associated metadata
void archive::explore(char * header, float ** metadata, profile * std = NULL){
  archive * original = new archive(this);
  fprintf(stderr,"Exploring %s with %s\n",this->filename,header);
  cpgeras();
  dialog * d = new dialog();
  int QUIT=d->addbutton(0.012,0.01,"QUIT");
  int PLOT=d->addbutton(0.012,0.06,"PLOT");
  int PURGE=d->addbutton(0.012,0.11,"PURGE");
  int MOVIE=d->addbutton(0.012,0.16,"MOVIE");
  int GRAY=d->addbutton(0.012,0.21,"GRAYSCALE");
  int FREQVSPHASE=d->addbutton(0.012,0.26,"FREQVSPHASE");
  int TIMEVSPHASE=d->addbutton(0.012,0.31,"TIMEVSPHASE");
  int MOWLAWN=d->addbutton(0.012,0.36,"MOW LAWN");
  int SEEKBOUNDS=d->addcheck(0.09,0.055,"SeekBounds?",1);
  int MAXHOLD=d->addcheck(0.09,0.1,"MaxHold",1);
  int SHOWZAPPED=d->addcheck(0.09,0.145,"ShowZeroWeights?",1);
  int MEDIANS=d->addcheck(0.09,0.19,"Plot medians",0);
  int PR0=d->addplotregion(0.3,0.99,0.1,0.6);
  int PR1=d->addplotregion(0.3,0.99,0.7,0.95);
  int PR3=d->addplotregion(0.1,0.25,0.7,0.95);
  float x,y;
  int plotno=0;
  int currentplottype=0;
  char ans;
  int button;
  float * chans = new float[nchan*nsubint];
  float * subs = new float[nchan*nsubint];
  float * weights = new float[nchan*nsubint];
  float * medians = new float[nchan];
  float * localmedians = new float[nchan];
  float * mads = new float[nchan];
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++)
      {
	chans[i*nchan+j]=(float)j;
	subs[i*nchan+j]=(float)i;
	weights[i*nchan+j]=this->subints[i]->polnprofiles[j]->profiles[0]->weight;
      }
    //fprintf(stderr,"calling parse\n");
  int rows = parse(header," ");
  float ** plotdata = new float * [rows+3];
  char ** titles = new char * [rows+3];
  char wordchannel[]="channel";
  char wordsubint[]="subint";
  char wordweight[]="weights";
  //printf("Copying %d data pointers\n",rows);
  for (int i=0;i<rows;i++)
    plotdata[i]=metadata[i];
  plotdata[rows]=chans;
  plotdata[rows+1]=subs;
  plotdata[rows+2]=weights;

  //printf("populating the titles\n");
  char * token = strtok(header," ");
  titles[0]=new char[strlen(token)+1];
  //printf("length of token %d %s\n",strlen(token),token);
  strcpy(titles[0],token);
  for (int i=1;i<rows;i++){
    //printf("Getting next token %d\n",i);
    token=strtok(NULL," ");
    //printf("length of token %d %s\n",strlen(token),token);
    titles[i]=new char[strlen(token)+1];
    strcpy(titles[i],token);
  }
  titles[rows]=wordchannel;
  titles[rows+1]=wordsubint;
  titles[rows+2]=wordweight;
  rows+=3;
  char space[]=" ";
  // create a row of radio controls
  for (int i=0;i<rows;i++){
    int value;
    if (i==0) value = 1; else value=0;
    d->addradio(0.07,0.5+(float)i/25.0,titles[i],value,0);
  }
  for (int i=0;i<rows;i++){
    int value;
    if (i==1) value = 1; else value=0;
    d->addradio(0.04,0.5+(float)i/25.0,space,value,1);
  }
  for (int i=0;i<rows;i++){
    int value;
    if (i==2) value = 1; else value=0;
    d->addradio(0.01,0.5+(float)i/25.0,space,value,2);
  }
  d->addstaticText(0.005, 0.5+(float)rows/25.0,"x",1,1.5,0.0);
  d->addstaticText(0.035, 0.5+(float)rows/25.0,"y",1,1.5,0.0);
  d->addstaticText(0.065, 0.5+(float)rows/25.0,"z",1,1.5,0.0);
  d->draw();

  // compute the ideal bandpass and pronouce
  float * bp = new float[nchan];
  fixbandpass(plotdata[rows-4],nchan,nsubint,11,1.1,0.25,bp);

  while ((button=d->manage(&x,&y,&ans,&plotno))!=QUIT){
    if (ans=='Q') break;
    int rx=d->groupon(2);
    int ry=d->groupon(1);
    int rz=d->groupon(0);
    int oldrx,oldry,oldrz;
    if (oldrx!=rx || oldry!=ry || oldrz!=rz){
      oldrx=rx; oldry=ry; oldrz=rz;
      d->checks[SEEKBOUNDS].on=1;
      d->update();
    }
    //printf("Seek Bounds is %d\n",d->checks[SEEKBOUNDS].on);
    
    if (plotno==0 && ans=='x'){
      printf("x hit x=%f y=%f plotnumber = %d\n",x,y,plotno);
      d->plotregions[plotno].reset();
      cpgsci(8);
      cpgmove(x,d->plotregions[plotno].yminworld);
      cpgdraw(x,d->plotregions[plotno].ymaxworld);
      cpgsci(1);
      float x1=x;
      cpgcurs(&x,&y,&ans);
      float x2=x;
      if (x2<x1){
	float temp=x1; x1=x2; x2=temp;
      }
      d->plotregions[0].set(x1,x2,d->plotregions[plotno].yminworld,
			    d->plotregions[plotno].ymaxworld);
      d->checks[SEEKBOUNDS].on=0;
      button=currentplottype;
    }
    if (plotno==0 && ans=='z'){
      printf("z hit x=%f y=%f plotnumber = %d\n",x,y,plotno);
      d->plotregions[plotno].reset();
      cpgsci(8);
      float x1=x,y1=y,x2,y2;
      char ans2;
      cpgband(2,1,x1,y1,&x2,&y2,&ans2);
      if (y2<y1){
	float temp=y2; y2=y1; y1=temp;
      }
      if (x2<x1){
	float temp=x2; x2=x1; x1=temp;
      }
      d->plotregions[0].set(x1,x2,y1,y2);
      d->checks[SEEKBOUNDS].on=0;
      button=currentplottype;
    }
    if (plotno==0 && ans=='u'){
      printf("z hit x=%f y=%f plotnumber = %d\n",x,y,plotno);
      d->plotregions[plotno].reset();
      d->checks[SEEKBOUNDS].on=1;
      d->update();
      button=currentplottype;
    }
    if (plotno==0 && ans=='y'){
      printf("y hit x=%f y=%f plotnumber = %d\n",x,y,plotno);
      d->plotregions[plotno].reset();
      cpgsci(8);
      cpgmove(d->plotregions[plotno].xminworld,y);
      cpgdraw(d->plotregions[plotno].xmaxworld,y);
      cpgsci(1);
      float y1=y;
      cpgcurs(&x,&y,&ans);
      float y2=y;
      if (y2<y1){
	float temp=y1; y1=y2; y2=temp;
      }
      d->plotregions[0].set(d->plotregions[plotno].xminworld,
		     d->plotregions[plotno].xmaxworld,y1,y2);
      d->checks[SEEKBOUNDS].on=0;
      button=currentplottype;
    }
    if (plotno==0 && ans=='q'){
      //printf("x=%f y=%f plotnumber = %d\n",x,y,plotno);
      d->plotregions[0].reset();
      float dx = d->plotregions[plotno].xmaxworld-d->plotregions[plotno].xminworld;
      float dy = d->plotregions[plotno].ymaxworld-d->plotregions[plotno].yminworld;
      //printf("dx=%f dy=%f\n",dx,dy);
      int ptofinterest = closest(nsubint*nchan,x,y,plotdata[rx],plotdata[ry],
				 dx,dy);
      d->plotregions[0].reset();
      cpgsci(2);
      cpgpt1(x,y,2);
      cpgsci(6);
      cpgpt1(plotdata[rx][ptofinterest],plotdata[ry][ptofinterest],17);
      cpgsci(1);
      d->plotregions[1].erase();
      d->plotregions[1].reset();
      if (std!=NULL)
      subints[(int)subs[ptofinterest]]->
	polnprofiles[(int)chans[ptofinterest]]->profiles[0]->plotresidual(std);
      else
	subints[(int)subs[ptofinterest]]->polnprofiles[(int)chans[ptofinterest]]->profiles[0]->plot(1);
      printf("Report for point %d\n",ptofinterest);
      printf("-----------------------------------------\n");
      for (int i=0;i<rows;i++){
	printf("%d %s %f\n",i,titles[i],plotdata[i][ptofinterest]);
      }
    }
    if (button==QUIT) printf("QUIT\n");
    if (button==PURGE){
      printf("PURGE\n");
      // if a point is in the plot set the weight to zero
      int count=0;
      for (size_t i=0;i<nsubint*nchan;i++){
	float xcoord=plotdata[rx][i];
	float ycoord=plotdata[ry][i];
	if (d->plotregions[0].insideworld(xcoord,ycoord)){
	    weights[i]=0.0;
	    size_t si=i/nchan;
	    size_t ci=i%nchan;
	    this->subints[si]->polnprofiles[ci]->set_weight(0.0);
	    //printf("Deleting subint %zu channel %zu\n",si,ci);
	    count++;
	  }
      }
      printf("Purged %d points\n",count);
      button=currentplottype;
      d->checks[SEEKBOUNDS].on=1;
    }
    if (button==MOWLAWN){
      printf("MOWLAWN\n");
      d->plotregions[0].reset();
      // if a point is more than 10% above the running median zap
      int count=0;
      cpgsci(0);
      for (size_t i=0;i<nsubint*nchan;i++){
	float xcoord=plotdata[rx][i];
	float ycoord=plotdata[ry][i];
	if (ycoord>1.1*bp[i%nchan]){
	    cpgpt1(xcoord,ycoord,1);
	    weights[i]=0.0;
	    size_t si=i/nchan;
	    size_t ci=i%nchan;
	    this->subints[si]->polnprofiles[ci]->set_weight(0.0);
	    //printf("Deleting subint %zu channel %zu\n",si,ci);
	    count++;
	  }
      }
      printf("Purged %d points\n",count);
      usleep(3000000);
      button=currentplottype;
      d->checks[SEEKBOUNDS].on=1;
      cpgsci(1);
    }
    if (button==PLOT) {
      currentplottype=PLOT;
      printf("PLOT\n");
      d->plotregions[0].erase();
      //d->plotregions[0].set(0,(float)nchan,0,(float)nsubint);
      //printf("%d %d %d\n",rx,ry,rz);
      // seek the bounds
      if (d->checks[SEEKBOUNDS].on) {
        float xmin5 = themin5(nchan*nsubint,plotdata[rx]);
        float xmax5 = themax5(nchan*nsubint,plotdata[rx]);
        float ymin5 = themin5(nchan*nsubint,plotdata[ry]);
        float ymax5 = themax5(nchan*nsubint,plotdata[ry]);
        //printf("%f %f %f %f\n",xmin5,xmax5,ymin5,ymax5);
        d->plotregions[0].set(xmin5,xmax5,ymin5,ymax5);
      }
      d->plotregions[0].reset();
      cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
      cpglab(titles[rx],titles[ry]," ");
      cpgpt(nchan*nsubint,plotdata[rx],plotdata[ry],1);
      d->checks[SEEKBOUNDS].on=0;
      float xdat[nchan];
      if (d->checks[MEDIANS].on){
	printf("Plotting medians and local medians\n");
	cpgsci(2);
	for (size_t i=0;i<nchan;i++){
	  xdat[i]=float(i);
	  float ypt=medianstride(plotdata[ry],nchan,nsubint,i);
	  float ypt2=madstride(plotdata[ry],nchan,nsubint,i);
	  medians[i]=ypt;
	  mads[i]=ypt2;
	  cpgpt1(xdat[i],ypt,17);
	  cpgmove(xdat[i],ypt-ypt2);
	  cpgdraw(xdat[i],ypt+ypt2);
	}
	local_medians((int)nchan,21,medians,localmedians);
	cpgsci(8);
	cpgline((int)nchan,xdat,localmedians);
	cpgsci(6);
	cpgline((int)nchan,xdat,bp);
	cpgsci(1);
      }
    }
    if (button==MOVIE) {
      currentplottype=MOVIE;
      printf("MOVIE seekbounds is %d\n",d->checks[SEEKBOUNDS].on);
      //d->plotregions[0].set(0,(float)nchan,0,(float)nsubint);
      //printf("%d %d %d\n",rx,ry,rz);
      // seek the bounds
      d->plotregions[0].erase();
      if (d->checks[SEEKBOUNDS].on) {
        float xmin5 = themin5(nchan*nsubint,plotdata[rx]);
        float xmax5 = themax5(nchan*nsubint,plotdata[rx]);
        float ymin5 = themin5(nchan*nsubint,plotdata[ry]);
        float ymax5 = themax5(nchan*nsubint,plotdata[ry]);
        d->plotregions[0].set(xmin5,xmax5,ymin5,ymax5);
	d->checks[SEEKBOUNDS].on=0;
      }
      d->plotregions[0].reset();
      cpgsci(1);
      cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
      cpglab(titles[rx],titles[ry]," ");
      for (size_t i=0;i<nsubint;i++){
	//cpgpt(nchan,&plotdata[rx][i*nchan],&plotdata[ry][i*nchan],1);
	cpgbbuf();
	if (i!=0){  //erase the old line
	  cpgscr(0,0.0,0.0,0.0);
	  cpgsci(0);
	  cpgslw(7);
	  cpgline(nchan,&plotdata[rx][(i-1)*nchan],&plotdata[ry][(i-1)*nchan]);
	  cpgslw(1);
	}
	if (i!=0 && d->checks[MAXHOLD].on) { // show the maxhold
	  cpgscr(15,0.2,0.2,0.6);
	  cpgsci(15);
	  cpgslw(3);
	  cpgline(nchan,&plotdata[rx][(i-1)*nchan],&plotdata[ry][(i-1)*nchan]);
	  cpgsci(1);
	}
	cpgsci(3);
	cpgslw(7);
	cpgline(nchan,&plotdata[rx][i*nchan],&plotdata[ry][i*nchan]);
	cpgslw(1);
	usleep(200000);
	cpgebuf();
      }
    }
    if (button==GRAY) {
      printf("GRAY\n");
      currentplottype=GRAY;
      d->plotregions[0].erase();
      // seek the bounds
      if (d->checks[SEEKBOUNDS].on) {
        float xmin = 0.0;
        float xmax = float(nchan);
        float ymin = 0.0;
        float ymax = float(nsubint);
        d->plotregions[0].set(xmin,xmax,ymin,ymax);
      }
      d->plotregions[0].reset();
      cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
      cpglab("Channel","Subint","Grayscale");
      float tr[6];
      tr[0]=0.0-0.5/float(nchan);
      tr[1]=1.0; tr[2]=0.0; tr[3]=0.0-0.5; tr[4]=0.0; tr[5]=1.0;
      heat();
      cpgimag(plotdata[rz],nchan,nsubint,1,nchan,1,nsubint,
	      themin(nsubint*nchan,plotdata[rz]),
	      themax(nsubint*nchan,plotdata[rz]),
	      tr);
      d->checks[SEEKBOUNDS].on=0;
    }
    d->update();  // update the dialog
  }
  printf("leaving Explore with button %d pressed\n",button);
}

void archive::set_DM(double input_dm){
  dm=input_dm;
}

double archive::dedisperse(){
  return(this->dedisperse(dm));
}

double archive::dedisperse(double input_dm){
  double f = this->get_weighted_frequency();
  for (size_t i=0;i<nsubint;i++){
    double p = subints[i]->folding_period;
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->dedisperse(f,input_dm,p);
  }
  return(f);
}

double archive::dedisperse(double input_dm, double refFreq){
  double f = refFreq;
  for (size_t i=0;i<nsubint;i++){
    double p = subints[i]->folding_period;
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->dedisperse(f,input_dm,p);
  }
  return(f);
}

void archive::remove_baseline(size_t ipol, int nbaseline){
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->profiles[ipol]->remove_baseline(nbaseline);
}

// rotates by turns
void archive::rotateleft(float shift){
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->rotateleft(shift);
}

// rotates by bins
void archive::rotateleft(int ishift){
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->rotateleft(ishift);
}

void archive::dedrift(double dTurnsdDay, long double tref){
  if (dTurnsdDay==0.0) return;
  double shift;
#pragma omp parallel for private (shift)
  for (size_t i=0;i<nsubint;i++){
    shift = dTurnsdDay*(subints[i]->mjd-tref);
    for (size_t j=0;j<nchan;j++)
      subints[i]->polnprofiles[j]->rotateleft(float(shift));
  }
}

float archive::max(size_t ipol){
  int found_one=0;
  float the_max;
  for (size_t i=0;i<nsubint;i++){
    if (subints[i]->nonzapped()){
      float test=subints[i]->max(ipol);
      if (found_one==0){
	found_one=1;
	the_max=test;
      }
      if (test>the_max) the_max=test;
    }
  }
  if (found_one) return(the_max);
  else return(0.0);
}

float archive::min(size_t ipol){
  int found_one=0;
  float the_min;
  for (size_t i=0;i<nsubint;i++){
    if (subints[i]->nonzapped()){
      float test=subints[i]->min(ipol);
      if (found_one==0){
	found_one=1;
	the_min=test;
      }
      if (test<the_min) the_min=test;
    }
  }
  if (found_one) return(the_min);
  else return(0.0);
}

double archive::get_weighted_frequency(size_t chan){
  return(this->get_weighted_frequency(chan,0,nsubint));
}

double archive::get_weighted_frequency(size_t chan, size_t startsub, size_t nadd){
  float sum_weight=0.0;
  double sum_weight_freq=0.0;
  for (size_t i=startsub;i<startsub+nadd && i<nsubint;i++){
    float w = subints[i]->polnprofiles[chan]->profiles[0]->weight;
    double f = subints[i]->polnprofiles[chan]->profiles[0]->freq;
    sum_weight_freq+= f*w;
    sum_weight+=w;
  }
  if (sum_weight==0.0) return(0.0);
  else 
    return(sum_weight_freq/sum_weight);
}

double archive::get_weighted_frequency(){
  float sum_weight=0.0;
  double sum_weight_freq=0.0;
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++){
      float w = subints[i]->polnprofiles[j]->profiles[0]->weight;
      double f = subints[i]->polnprofiles[j]->profiles[0]->freq;
      sum_weight_freq+= f*w;
      sum_weight+=w;
    }
  if (sum_weight==0.0) return(0.0);
  else 
    return(sum_weight_freq/sum_weight);
}

void archive::fadd(int ntoadd){
  int nadd=ntoadd;
  if (ntoadd==0) {
    nadd=nchan;
  }
  if (nchan%nadd==0){
    for (int i=0;i<(int)nsubint;i++)
      subints[i]->fadd(nadd,dm);
    nchan=nchan/nadd;
  } else {
    fprintf(stderr,"archive::fadd %d ntoadd must divide nchan %d",ntoadd,(int)nchan);
  }
}

// assumes a pgplot window is open already
void archive::multiplot(){
  for (int i=0;i<(int)nsubint;i++)
    for (int j=0;j<(int)nchan;j++)
      {
	// set the plotting window viewport
	float dx = 1.0 / (float) nchan;
	float dy = 1.0 / (float) nsubint;
	float xmin = (float) j * dx;
	float xmax = xmin+dx;
        float ymin = (float) i * dy;
        float ymax = ymin+dy;
	cpgsvp(xmin,xmax,ymin,ymax);
	subints[i]->polnprofiles[j]->profiles[0]->plot(0);
      }
}

// assumes a pgplot window is open already
void archive::Splot(size_t isub, size_t ichan, int click){
  float dy = 1.0 / (float) npol * 0.75;
  for (size_t k=0;k<npol;k++){
	float xmin = 0.15;
	float xmax = 0.9;
        float ymin = (float) k * dy + 0.15;
        float ymax = ymin+dy;
	//cpgsvp(xmin,xmax,ymin,ymax);
	subints[isub]->polnprofiles[ichan]->profiles[k]->plot(1);
  }
}

// assumes a pgplot window is open already
void archive::Splot(int click){
  this->Splot(0,0,click);
}

// assumes a pgplot window is open already
void archive::Gplot(size_t ipol, int click){
  this->remove_baseline(ipol,nbin/8);
  float zmin=this->min(ipol);
  float zmax=this->max(ipol);
  float xmin=0.0;
  float xmax=1.0;
  float ymin=0.0;
  float ymax=(float)nchan;
  float wmax = subints[0]->get_max_weight();
  float * d = new float[nbin*nchan];
  if (d==NULL){
    fprintf(stderr,"archive::Gplot error allocating plot region for %zu bytes\n",
	    nbin*nchan*sizeof(float));
    exit(-1);
  }
  // fill with the minimum zvalue
  for (size_t isub=0;isub<nsubint;isub++){
    for (size_t i=0;i<nchan;i++)
      for (size_t j=0;j<nbin;j++)
	d[i*nbin+j]=0.0;
    for (size_t i=0;i<nchan;i++){
      if (subints[isub]->polnprofiles[i]->profiles[ipol]->weight!=0.0){
      //fprintf(stderr,"Copying subint %d of %d\n",(int) i, (int)nsubint);
	memcpy(&d[i*nbin],
	       subints[isub]->polnprofiles[i]->profiles[ipol]->amps,
	       nbin*sizeof(float));
	for (size_t ibin=0;ibin<nbin;ibin++)
	d[i*nbin+ibin]*=subints[isub]->polnprofiles[i]->profiles[ipol]->weight/wmax;
      }
    }
  float tr[6];
  tr[0]=0.0-0.5/float(nbin);
  tr[1]=1.0/float(nbin); tr[2]=0.0; tr[3]=0.0-0.5; tr[4]=0.0; tr[5]=1.0;
  //  cpgsvp(0.1,0.9,0.1,0.9);
  cpgswin(xmin,xmax,ymin,ymax);
  heat();
  cpgimag(d,(int)nbin,(int)nchan,1,(int)nbin,1,(int)nchan,
	  zmin,zmax,tr);
  cpgbox("BCST",0.0,0,"BCST",0.0,0);
  cpgmtxt("L",1.0,0.5,0.5,"Frequency Channel");
  //cpgmtxt("T",1.0,0.5,0.5,this->filename);
  }
  float x,y;
  char ans;
  if (click) cpgcurs(&x,&y,&ans);
  if (ans=='X') printf("%s %c\n",this->filename,ans);
}

// assumes a pgplot window is open already
void archive::Yplot(size_t ipol, int click){
  // scale data by the weight of each subint
  archive * copy = new archive(this);
  for (size_t i=0;i<nsubint;i++)
    copy->subints[i]->polnprofiles[0]->profiles[ipol]->scale(copy->subints[i]->polnprofiles[0]->profiles[ipol]->weight);
  copy->remove_baseline(ipol,nbin/8);
  float zmin=copy->min(ipol);
  float zmax=copy->max(ipol);
  float xmin=0.0;
  float xmax=1.0;
  float ymin=0.0;
  float ymax=(float)nsubint-1.0;
  float * d = new float[copy->nbin*copy->nsubint];
  if (d==NULL){
    fprintf(stderr,"archive::Yplot error allocating plot region for %zu bytes\n",
	    nbin*nsubint*sizeof(float));
    exit(-1);
  }
  // fill with the minimum zvalue
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nbin;j++)
      d[i*nbin+j]=zmin;
  for (size_t i=0;i<nsubint;i++){
    if (copy->subints[i]->nonzapped()){
      //fprintf(stderr,"Copying subint %d of %d\n",(int) i, (int)nsubint);
      memcpy(&d[i*nbin],
	     copy->subints[i]->polnprofiles[0]->profiles[ipol]->amps,
	     copy->nbin*sizeof(float));
    }
  }
  float tr[6];
  tr[0]=0.0-0.5/float(nbin);
  tr[1]=1.0/float(nbin); tr[2]=0.0; tr[3]=0.0-0.5; tr[4]=0.0; tr[5]=1.0;
  cpgswin(xmin,xmax,ymin,ymax);
  heat();
  cpgimag(d,(int)nbin,(int)nsubint,1,(int)nbin,1,(int)nsubint,
	  zmin,zmax,tr);
  cpgbox("BCST",0.0,0,"BCST",0.0,0);
  cpgmtxt("L",1.0,0.5,0.5,"Time index");
}

int archive::showHeader(int verbose=0){
  if (verbose){
    printf("Header for file: %s\n",filename);
    printf("           nbin: %d\n",(int)nbin);
    printf("           npol: %d\n",(int)npol);
    printf("          nchan: %d\n",(int)nchan);
    printf("        nsubint: %d\n",(int)nsubint);
    printf("      Frequency: %10.5lf\n",freq);
    printf(" bandwidth(MHz): %f\n",bw);
    printf("             DM: %lf\n",dm);
    printf("         Source: %-15.15s\n",source);
    printf("      Telescope: %-15.15s\n",telescope);
    printf("       RA (rad): %lf\n",ra);
    printf("      Dec (rad): %lf\n",dec);
  } else {
    printf("%s ",filename);
    printf("%d ",(int)nbin);
    printf("%d ",(int)npol);
    printf("%d ",(int)nchan);
    printf("%d ",(int)nsubint);
    printf("%10.5lf ",freq);
    printf("%f ",bw);
    printf("%lf ",dm);
    printf("%-15.15s ",source);
    printf("%-15.15s ",telescope);
    printf("%lf ",ra);
    printf("%lf \n",dec);
  }
  return(0);
}

int archive::showIntHeader(){
  for (size_t i=0;i<nsubint;i++) {
    printf("%zu ",i);
    printf("%s ",filename);
    printf("%20.14Lf ",subints[i]->mjd);
    printf("%10.6lf ",subints[i]->duration);
    printf("%14.12lf ",subints[i]->folding_period);
    printf("%8.4f ",subints[i]->az*180.0/M_PI);
    printf("%8.4f ",subints[i]->zen*180.0/M_PI);
    printf("%8.4f ",subints[i]->parallactic*180.0/M_PI);
    printf("%8.4f\n",subints[i]->meridianDeg());
  }
  return(0);
}

int archive::showWeightHeader(){
  for (size_t i=0;i<nsubint;i++)
    for (size_t j=0;j<nchan;j++)
      fprintf(stderr,"isub: %d ichan: %d weight: %f freq: %lf\n",(int)i,(int)j,
	      subints[i]->polnprofiles[j]->profiles[0]->weight,
	      subints[i]->polnprofiles[j]->profiles[0]->freq);
  return(0);
}

// add together the subints in groups of nadd
// if nadd==0 then add them all!
// requires dedispersion to a common frequency per polnprofile
// before adding
void archive::tscrunch(size_t nadd){

  //fprintf(stderr,"tscrunching by a factor %zu (0=whole thing)\n",nadd);
  // special case
  // still need to update mjd, az, ze, para etc.
  if (nadd==0){
    for (size_t i=1; i<nsubint;i++)
      subints[0]->duration+=subints[i]->duration;
    for (size_t j=0; j<nchan;j++) {
      // get the weighted freq for each channel
      double weighted_freq=this->get_weighted_frequency(j);
      // then dedisperse them
      double fp;
      #pragma omp parallel for private(fp)
      for (size_t i=0; i<nsubint;i++){
	double fp=subints[i]->folding_period;
	subints[i]->polnprofiles[j]->dedisperse(weighted_freq,dm,fp);
      }
      for (size_t i=1; i<nsubint;i++){
      #pragma omp parallel for 
	for (size_t k=0;k<npol;k++)
	  subints[0]->polnprofiles[j]->profiles[k]->add(
			      subints[i]->polnprofiles[j]->profiles[k]);
      }
    }
    nsubint = 1;
  }
    else {
  // needs the durations computed plus the appropriate centre frequencies
    size_t ntargets = nsubint/nadd;
    if ( nsubint - (ntargets*nadd) != 0 ) ntargets++;
    //fprintf(stderr,"tscrunching to %zu targets\n",ntargets);
    // first determine the new frequencies of all the subints for each channel
    double cfs[nchan][ntargets];
    for (size_t j=0; j<nchan;j++)
      for (size_t itarget=0; itarget<ntargets;itarget++){
	//fprintf(stderr,"Getting weighted_freq \n");
	  cfs[j][itarget] = get_weighted_frequency(j,itarget*nadd,nadd);
      }
    // now dedisperse them all in the most massively parallel way
    #pragma omp parallel for
    for (size_t i=0; i<nsubint;i++)
      for (size_t j=0;j<nchan;j++){
	//fprintf(stderr,"dedispersing \n");
	subints[i]->polnprofiles[j]->dedisperse(
			   cfs[j][i/nadd],dm,subints[i]->folding_period);
      }
    // now add them all up - being careful not to overwrite anything.
    #pragma omp parallel for
    for (size_t j=0;j<nchan;j++)
      #pragma omp parallel for
      for (size_t k=0;k<npol;k++)
	for (size_t itarget=0;itarget<ntargets;itarget++){
	  //fprintf(stderr,"target %zu of %zu total sub %zu\n",itarget,ntargets,nsubint);
	  subints[itarget]->polnprofiles[j]->profiles[k]->copy(
	       subints[itarget*nadd]->polnprofiles[j]->profiles[k]);
	  for (size_t i=1;i<nadd && i+itarget*nadd<nsubint;i++){
	    subints[itarget]->polnprofiles[j]->profiles[k]->add(
	       subints[itarget*nadd+i]->polnprofiles[j]->profiles[k]);
	  }
	}
    // Now update the tricky stuff in a precise way
    // could probably be parallelised too
    for (size_t itarget=0;itarget<ntargets;itarget++){
	  long double mjd_new = get_median_mjd(itarget*nadd,nadd);
	  double duration_new = get_total_duration(itarget*nadd,nadd);
	  double fp_new = get_folding_period_median(itarget*nadd,nadd);
	  subints[itarget]->mjd=mjd_new;
	  subints[itarget]->duration=duration_new;
	  subints[itarget]->folding_period=fp_new;
	  // need to add in weighted az, zen, para etc too
    }
  nsubint = ntargets;
    }
}

void archive::pscrunch(int plan){
  if (npol!=1){
    #pragma omp parallel for
    for (size_t i=0; i<nsubint;i++)
      subints[i]->pscrunch(plan);
    npol=1;
  }
}

void archive::bscrunch(int nadd){
    #pragma omp parallel for
    for (size_t i=0; i<nsubint;i++)
      subints[i]->bscrunch(nadd);
    nbin=nbin/nadd;
}

#define HDRSIZE 80
#define INTHDRINCREMENT 44 
#define CHANHDRINCREMENT 12

// load archive from filename
// populate data and handy classes
archive::archive(char * fname){
  struct stat buffer;
  FILE * fptr = fopen(fname,"r");
  if (fptr==NULL){
    fprintf(stderr,"archive::archive Error opening file %s\n",fname);
    exit(-1);
  }
  fstat(fileno(fptr),&buffer);
  //fprintf(stderr,"archive::archive Opened file %s that has %zu bytes in it\n",fname,buffer.st_size);
  filename = new char [strlen(fname)+1];
  strcpy(filename,fname);
  unsigned char header[HDRSIZE];
  fread(header,sizeof(unsigned char),HDRSIZE,fptr);

  // nbin, npol, nchan, nsubint
  int testint;
  memcpy(&testint,&header[0],sizeof(int));
  nbin = (size_t) testint;
  memcpy(&testint,&header[4],sizeof(int));
  npol = (size_t) testint;
  memcpy(&testint,&header[8],sizeof(int));
  nchan = (size_t) testint;
  memcpy(&testint,&header[12],sizeof(int));
  nsubint = (size_t) testint;

  // sanity checks

  if (nbin<1 || nbin>16384 || nchan<1 || nchan>4096 || nsubint<1 || npol<1 || npol>4){
    fprintf(stderr,"This file (%s) appears not to be a sensible raw observation file\n",
	    fname);
    fprintf(stderr,"nbin %zu nchan %zu npol %zu nsubint %zu\n",nbin,nchan,npol,nsubint);
    exit(-1);
  }

  original_nbin = nbin;
  original_npol = npol;
  original_nchan = nchan;
  original_nsubint = nsubint;

  // freq, bw, dm
  memcpy(&freq,&header[16],sizeof(float));
  memcpy(&bw,&header[20],sizeof(float));
  memcpy(&dm,&header[24],sizeof(double));

  // The source and telescope name
  char * token = NULL;
  token = strtok((char *) &header[32]," ");  
  strncpy(source,token,16);
  token = strtok((char *) &header[48]," ");  
  strncpy(telescope,token,16);

  // ra, dec
  memcpy(&ra,&header[64],sizeof(double));
  memcpy(&dec,&header[72],sizeof(double));

  //fprintf(stderr,"Successfully read in the header of length %zu bytes\n",(size_t)HDRSIZE);

  // Now read in the integration header
  unsigned char * intbuffer = new unsigned char [(size_t)INTHDRINCREMENT * nsubint];
  fread(intbuffer,sizeof(unsigned char),(size_t) INTHDRINCREMENT * nsubint,fptr);
  //fprintf(stderr,"Successfully read in the subint headers total %zu\n",INTHDRINCREMENT * nsubint);

  // Now read in the profile weights and frequencies
  unsigned char * chanbuffer = new unsigned char [(size_t)CHANHDRINCREMENT * nsubint * nchan];
  fread(chanbuffer,sizeof(unsigned char),(size_t)CHANHDRINCREMENT * nsubint * nchan,fptr);
  //  fprintf(stderr,"Successfully read in the profile weights and freqs total %zu\n",
  //  (size_t) CHANHDRINCREMENT * nsubint * nchan);

  data = new float[(size_t)nbin*(size_t)npol*(size_t)nchan*(size_t)nsubint];
  if (data==NULL){
    fprintf(stderr,"archive::archive Error creating array data %zu %zu %zu %zu\n",nbin,
	    npol,nchan,nsubint);
  }
  //fprintf(stderr,"archive::archive created RAM for archive %s\n",filename);
  fread(data,sizeof(float),(size_t)nbin*(size_t)npol*(size_t)nchan*(size_t)nsubint,fptr);
  //fprintf(stderr,"archive::archive filled RAM for archive %s\n",filename);
  fclose(fptr);
  //fprintf(stderr,"archive::archive closed file %s\n",filename);
  // create and populate all of the pointers for easy use here
  size_t dprof = nbin;
  size_t dpolnprof = npol * dprof;
  size_t dsubint = dpolnprof * nchan;
  subints = new subint * [nsubint];     // an array of pointers to the subint pointers   
  for (size_t i=0;i<nsubint;i++){
    subints[i] = new subint(nchan,npol);
    subints[i]->populate((float*)&intbuffer[i*(size_t)INTHDRINCREMENT]);
    for (size_t j=0;j<nchan;j++){
      subints[i]->polnprofiles[j] = new polnprofile(nbin,npol);
      for (size_t k=0;k<npol;k++){
	subints[i]->polnprofiles[j]->profiles[k] = new profile(nbin,&data[i*dsubint+j*dpolnprof+k*dprof]);
      }
      float * w = (float *) &chanbuffer[sizeof(float)*(i*nchan+j)];
      double * f = (double *) &chanbuffer[sizeof(float)*nsubint*nchan+(i*nchan+j)*sizeof(double)];
      subints[i]->polnprofiles[j]->set_weight(w); 
      subints[i]->polnprofiles[j]->set_freq(f);
    }
  }
//  fprintf(stderr,"successfully unpacked %s\n",filename);
}

archive::archive(char * abuffer, size_t bufferlen){
  filename = NULL;
  char * header = abuffer;
  // nbin, npol, nchan, nsubint
  int testint;
  memcpy(&testint,&header[0],sizeof(int));
  nbin = (size_t) testint;
  memcpy(&testint,&header[4],sizeof(int));
  npol = (size_t) testint;
  memcpy(&testint,&header[8],sizeof(int));
  nchan = (size_t) testint;
  memcpy(&testint,&header[12],sizeof(int));
  nsubint = (size_t) testint;

  original_nbin = nbin;
  original_npol = npol;
  original_nchan = nchan;
  original_nsubint = nsubint;

  // freq, bw, dm
  memcpy(&freq,&header[16],sizeof(float));
  memcpy(&bw,&header[20],sizeof(float));
  memcpy(&dm,&header[24],sizeof(double));

  // The source and telescope name
  strncpy(source,(char *)&header[32],16);
  strncpy(telescope,(char *)&header[48],16);

  // ra, dec
  memcpy(&ra,&header[64],sizeof(double));
  memcpy(&dec,&header[72],sizeof(double));

  //fprintf(stderr,"Successfully read in the header of length %zu bytes\n",(size_t)HDRSIZE);
  // Now read in the integration header
  char * intbuffer = (char *) &abuffer[(size_t)HDRSIZE];

  // Now read in the profile weights and frequencies
  char * chanbuffer = (char *) &abuffer[(size_t)HDRSIZE+
							(size_t)INTHDRINCREMENT * nsubint];
  data = (float *) &abuffer[(size_t)HDRSIZE+(size_t)INTHDRINCREMENT*nsubint+
			    (size_t)CHANHDRINCREMENT*nsubint*nchan];

  // create and populate all of the pointers for easy use here
  size_t dprof = nbin;
  size_t dpolnprof = npol * dprof;
  size_t dsubint = dpolnprof * nchan;
  subints = new subint * [nsubint];     // an array of pointers to the subint pointers   
  for (size_t i=0;i<nsubint;i++){
    subints[i] = new subint(nchan,npol);
    subints[i]->populate((float*)&intbuffer[i*(size_t)INTHDRINCREMENT]);
    for (size_t j=0;j<nchan;j++){
      subints[i]->polnprofiles[j] = new polnprofile(nbin,npol);
      for (size_t k=0;k<npol;k++){
	subints[i]->polnprofiles[j]->profiles[k] = new profile(nbin,&data[i*dsubint+j*dpolnprof+k*dprof]);
      }
      float * w = (float *) &chanbuffer[sizeof(float)*(i*nchan+j)];
      double * f = (double *) &chanbuffer[sizeof(float)*nsubint*nchan+(i*nchan+j)*sizeof(double)];
      subints[i]->polnprofiles[j]->set_weight(w); 
      subints[i]->polnprofiles[j]->set_freq(f);
    }
  }
  //fprintf(stderr,"successfully loaded %s\n",filename);
}

archive::archive(archive * from){
  filename = new char[strlen(from->filename)+1];
  strcpy(filename,from->filename);
  nbin = from->nbin;
  npol = from->npol;
  nchan = from->nchan;
  nsubint = from->nsubint;

  original_nbin = nbin;
  original_npol = npol;
  original_nchan = nchan;
  original_nsubint = nsubint;

  // freq, bw, dm

  freq=from->freq;
  bw=from->bw;
  dm=from->dm;

  // The source and telescope name
  strncpy(source,from->source,16);
  strncpy(telescope,from->telescope,16);

  // ra, dec
  ra=from->ra;
  dec=from->dec;
  //  fprintf(stderr,"Successfully read in the header of length %zu bytes\n",(size_t)HDRSIZE);

  data = new float[(size_t)nbin*(size_t)npol*(size_t)nchan*(size_t)nsubint];
  if (data==NULL){
    fprintf(stderr,"archive::archive Error creating array data %zu %zu %zu %zu\n",nbin,
	    npol,nchan,nsubint);
  }
  subints = new subint * [nsubint];     // an array of pointers to the subint pointers   
  for (size_t i=0;i<nsubint;i++){
    subints[i] = new subint(nchan, npol);
    subints[i]->copyhdr(from->subints[i]);
    for (size_t j=0;j<nchan;j++){
      subints[i]->polnprofiles[j] = new polnprofile(nbin,npol);
      for (size_t k=0;k<npol;k++){
	// copy the amps into data
	memcpy(&data[nbin*nchan*npol*i+nbin*npol*j+k*nbin],
	       &(from->subints[i]->polnprofiles[j]->profiles[k]->amps[0]),
	       sizeof(float)*nbin);
	subints[i]->polnprofiles[j]->profiles[k] = new profile(nbin,
	      &data[nbin*nchan*npol*i+nbin*npol*j+k*nbin]);
      }
      float w = from->subints[i]->polnprofiles[j]->profiles[0]->weight;
      double f = from->subints[i]->polnprofiles[j]->profiles[0]->freq;
      subints[i]->polnprofiles[j]->set_weight(&w); 
      subints[i]->polnprofiles[j]->set_freq(&f);
    }
  }
  //fprintf(stderr,"successfully copied archive\n");
}

void archive::refill(archive * from){
  strcpy(filename,from->filename);
  nbin = from->nbin;
  npol = from->npol;
  nchan = from->nchan;
  nsubint = from->nsubint;

  original_nbin = nbin;
  original_npol = npol;
  original_nchan = nchan;
  original_nsubint = nsubint;

  // freq, bw, dm

  freq=from->freq;
  bw=from->bw;
  dm=from->dm;

  // The source and telescope name
  strncpy(source,from->source,16);
  strncpy(telescope,from->telescope,16);

  // ra, dec
  ra=from->ra;
  dec=from->dec;
  //  fprintf(stderr,"Successfully read in the header of length %zu bytes\n",(size_t)HDRSIZE);
  for (size_t i=0;i<nsubint;i++){
    subints[i]->copyhdr(from->subints[i]);
    for (size_t j=0;j<nchan;j++){
      for (size_t k=0;k<npol;k++){
	// copy the amps into data
	memcpy(&(subints[i]->polnprofiles[j]->profiles[k]->amps[0]),
	       &(from->subints[i]->polnprofiles[j]->profiles[k]->amps[0]),
	       sizeof(float)*nbin);
      }
      float w = from->subints[i]->polnprofiles[j]->profiles[0]->weight;
      double f = from->subints[i]->polnprofiles[j]->profiles[0]->freq;
      subints[i]->polnprofiles[j]->set_weight(&w); 
      subints[i]->polnprofiles[j]->set_freq(&f);
    }
  }
  //fprintf(stderr,"successfully copied archive\n");
}

void archive::incoadd(archive * from){
  if ((nbin == from->nbin) &&
      (npol == from->npol) &&
      (nchan == from->nchan) &&
      (nsubint == from->nsubint) &&
      (freq==from->freq) &&
      (bw==from->bw) &&
      (dm==from->dm) &&
      (strncmp(source,from->source,16)==0) &&
      (strncmp(telescope,from->telescope,16)==0) &&
      (ra == from->ra) &&
      (dec == from->dec)){
    for (size_t i=0;i<nsubint;i++)
      for (size_t j=0;j<nchan;j++)
	for (size_t k=0;k<npol;k++)
	// copy the amps into data
	   subints[i]->polnprofiles[j]->profiles[k]->add(
	     from->subints[i]->polnprofiles[j]->profiles[k]);
    //fprintf(stderr,"successfully summed archive %s into %s\n",
    //	    from->filename,this->filename);
  } else{
    fprintf(stderr,"archives %s and %s are incompatible\n",from->filename,this->filename);
    from->showHeader();
    this->showHeader();
  }
}

int archive::unload(){
  return(this->unload(this->filename));
}

// needs testing
int archive::unload(char * root, char * extension){
  printf("archive::unload(root,ext) >%s<\n",root);
  printf("archive::unload(root,ext) >%s<\n",extension);

  char newname[strlen(root)+strlen(extension)+1];
  if (strlen(extension)==0){
    strcpy(newname,root);
  }else{
    sprintf(newname,"%s%s",root,extension);
    newname[strlen(root)+strlen(extension)]='\0';
  }
  return (this->unload(newname));
}

int archive::unload(char * newname){
  printf("archive::unload newname >%s<\n",newname);
  FILE * fptr = fopen(newname,"w");
  if (fptr==NULL){
    fprintf(stderr,"archive::unload Error opening file %s\n",newname);
    exit(-1);
  }
  unsigned char header[HDRSIZE];

  // nbin, npol, nchan, nsubint
  int testint;
  testint = int(nbin);
  memcpy(&header[0],&testint,sizeof(int));
  testint = int(npol);
  memcpy(&header[4],&testint,sizeof(int));
  testint = int(nchan);
  memcpy(&header[8],&testint,sizeof(int));
  testint = int(nsubint);
  memcpy(&header[12],&testint,sizeof(int));

  // freq, bw, dm
  memcpy(&header[16],&freq,sizeof(float));
  memcpy(&header[20],&bw,sizeof(float));
  memcpy(&header[24],&dm,sizeof(double));

  // Blank the write source and telescope name
  strncpy((char *) &header[32],"             ",16);
  strncpy((char *) &header[48],"             ",16);
  strncpy((char *) &header[32],source,16);
  strncpy((char *) &header[48],telescope,16);

  // ra, dec
  memcpy(&header[64],&ra,sizeof(double));
  memcpy(&header[72],&dec,sizeof(double));

  fwrite(header,sizeof(unsigned char),HDRSIZE,fptr);
  //  fprintf(stderr,"Successfully wrote out the header of length %zu bytes\n",(size_t)HDRSIZE);

  // Now create room for the integration headers
  unsigned char * intbuffer = new unsigned char [(size_t)INTHDRINCREMENT * nsubint];
  for (size_t i=0;i<nsubint;i++)
    subints[i]->depopulate((float*)&intbuffer[i*(size_t)INTHDRINCREMENT]);
  fwrite(intbuffer,sizeof(unsigned char),(size_t) INTHDRINCREMENT * nsubint,fptr);
  //fprintf(stderr,"Successfully wrote out the subint headers total %zu\n",INTHDRINCREMENT * nsubint);

  // Now retrieve then write out the profile weights and frequencies
  unsigned char * chanbuffer = new unsigned char [(size_t)CHANHDRINCREMENT * nsubint * nchan];
  if (chanbuffer==NULL){
    fprintf(stderr,"archive::unload Error: cannot allocate %zu bytes for weight-freq header\n",
	    (size_t)CHANHDRINCREMENT * nsubint * nchan);
    return(-1);
  }
  size_t dprof = nbin;
  size_t dpolnprof = npol * dprof;
  size_t dsubint = dpolnprof * nchan;
  for (size_t i=0;i<nsubint;i++){
    for (size_t j=0;j<nchan;j++){
      float w = subints[i]->polnprofiles[j]->profiles[0]->weight;
      double f = subints[i]->polnprofiles[j]->profiles[0]->freq;
      memcpy(&chanbuffer[sizeof(float)*(i*nchan+j)],&w,sizeof(float));
      memcpy(&chanbuffer[sizeof(float)*nsubint*nchan+(i*nchan+j)*sizeof(double)],&f,sizeof(double));
    }
  }
  fwrite(chanbuffer,sizeof(unsigned char),(size_t)CHANHDRINCREMENT * nsubint * nchan,fptr);
  //fprintf(stderr,"Successfully wrote out the profile weights and freqs total %zu\n",
  //	  (size_t) INTHDRINCREMENT * nsubint);

  // Could save a lot of time and effort if original dimensions = current ones or
  // you have only tscrunched
  // Alternative is to just do one subint at a time using lots of fwrites of reasonable
  // dimension 
  float * newdata = new float[(size_t)nbin*(size_t)npol*(size_t)nchan*(size_t)nsubint];
  if (newdata==NULL){
    fprintf(stderr,"archive::unload Error creating array data %zu x %zu x %zu x %zu x sizeof(float)\n",nbin,
	    npol,nchan,nsubint);
  }
  //fprintf(stderr,"archive::unload created RAM for new archive %s\n",newname);
  // Now populate newdata before writing it out.
  // massively parallel memory sloshing
  for (size_t i=0;i<nsubint;i++){
    for (size_t j=0;j<nchan;j++){
      for (size_t k=0;k<npol;k++){
	memcpy(&newdata[i*dsubint +j*dpolnprof +k*dprof],
	       subints[i]->polnprofiles[j]->profiles[k]->amps,sizeof(float)*nbin);
      }
    }
  }
  fwrite(newdata,sizeof(float),(size_t)nbin*(size_t)npol*(size_t)nchan*(size_t)nsubint,fptr);
  fprintf(stderr,"archive::unloaded created archive %s\n",newname);
  fclose(fptr);
  fprintf(stderr,"successfully unloaded %s\n",newname);
  return(0);
}

int main(int argc, char ** argv){
  extern char *optarg;
  extern int optind;
  int c;
  static char usage[] = "usage: %s -0 -4 -a -b factor -c -d -e ext -f factor -g pgplotdev -h -k -o -n newpsrname -r bandpassfactor -s std -t tscrunch_factor -v -w \"minbins,maxbins,factor\" -x snr_max -z -A newname -B -C \"Gain,Trec\"-D DM -E edge -F -G -H -I -K -L -M maxang -N \"nx,ny\" -S -R width(turns) -P \"nsub,nchan\" -T -U -X -Y -Z\n";
  char * filename;
  archive * grand = NULL;
  int verbose=0, help=0;
  int plot4 = 0;
  int zap=0;
  int addarchs=0;
  int plot=0;
  int multiplot=0;
  int centre=0;
  int tscrunch=0;
  int Fscrunch=0;
  int fscrunch=0;
  int pscrunch=0;
  int bscrunch=0;
  int Nplot=1;
  int nxplot=1;
  int click=0;
  int nyplot=1;
  int hasplotdev=0;
  int dedisperse=0;
  int Yplot=0;
  int Gplot=0;
  int Splot=0;
  int computeBaseline=0;
  int Splot_sub=0;
  int Splot_chan=0;
  int showHeader=0;
  int showIntHeader=0;
  int showWeightHeader=0;
  int Tscrunch=0;
  int bfactor=1;
  int extension=0;
  int newpsr=0;
  int newdm=0;
  int hasstd=0;
  int mkstd=0;
  int psrchiveZap=0;
  int ascii=0;
  int Xplore=0;
  int getflux=0;
  float worst_SNR=6.0;
  float minlength=0.0;
  double input_dm;
  char * grandfilename=NULL;
  char * plotdev = NULL;
  char * ext_chars = NULL; // a silly length for an extension
  char * newpsrname=NULL;
  char stdfile[1000]; // a stupidly long filename
  float basew=0.3; // default
  float wstart = 0.01; //default
  float wend = 0.5; //default
  float wfactor = 1.5; //default of 3/2
  float istart = 0.01; //default
  float iend = 0.7; //default
  float ifactor = 1.5; //default of 3/2
  float irms = 0.0; //default
  float rfi_factor = 1.20; //default
  int fscrunch_factor=1;
  int tscrunch_factor=1;
  int plotcounter=0;
  int grandcount=0;
  int seekflag=1;
  int optimise=0;
  int ReWeight=0;
  float timeturns=0.0;
  float dmturns=0.0;
  float Gain=0.0;  // guilty until proven to have a value
  float MeerKATGain=3.364; // K/Jy
  float Trec=0.0; // detonates if not parsed
  float TCMB=2.725; // Kelvin
  int Cal=0;
  int maxMeridian=0;
  float maxMeridianAngle=180.0;
  int edgeZap=0;
  float edgeFrac=0.0;
  int gettoas=0;
  float ibasewidth=0.3;
  float psrbasewidth=0.3;
  int lowsnr=0;
  int nscan=0;

  while ((c=getopt(argc, argv, "04A:BC:D:E:FGHIKLM:P:N:O:Q:R:STUWXYZab:cde:f:g:hi:kn:opr:s:t:u:vw:x:z"))!= -1){
    switch (c) {
    case '0':
      seekflag=0;
      break;
    case '4':
      plot4=1;
      plot=1;
      break;
    case 'k':
      click=1;
      break;
    case 'n':
      newpsr=1;
      newpsrname = new char[strlen(optarg)+1];
      nscan = sscanf(optarg,"%s",newpsrname);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as the new pulsar name\n",optarg);
	exit(-1);
      }
      break;
    case 'L':
      lowsnr=1;
      break;
    case 'R':
      nscan = sscanf(optarg,"%f",&ibasewidth);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as RFInterference baseline width\n",optarg);
	exit(-1);
      }
      break;
    case 'Q':
      nscan = sscanf(optarg,"%f",&psrbasewidth);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as Pulsar baseline width\n",optarg);
	exit(-1);
      }
      break;
    case 'O':
      optimise = 1;
      nscan = comma2floatscan(optarg,&timeturns,&dmturns);
      if (nscan!=2){
	fprintf(stderr,"Error parsing %s as time and dm turns\n",optarg);
	exit(-1);
      }
      break;
    case 'C':
      Cal = 1;
      nscan = comma2floatscan(optarg,&Gain,&Trec);
      if (nscan!=2){
	fprintf(stderr,"Error parsing %s as Gain (K/Jy) and Trec (K)\n",optarg);
	exit(-1);
      }
      break;
    case 'A':
      addarchs=1;
      grandfilename = new char[strlen(optarg)+1];
      nscan = sscanf(optarg,"%s",grandfilename);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as grandfilename\n",optarg);
	exit(-1);
      }
      break;
    case 'B':
      computeBaseline=1;
      break;
    case 'N':
      Nplot=1;
      nscan = comma2intscan(optarg,&nxplot,&nyplot);
      if (nscan!=2){
	fprintf(stderr,"Error parsing %s as NXplot and NYplot\n",optarg);
	exit(-1);
      }
      break;
    case 'a':
      ascii=1;
      break;
    case 'c':
      centre=1;
      break;
    case 'i':
      nscan = comma3floatscan(optarg,&istart,&iend,&ifactor);
      if (nscan!=3){
	fprintf(stderr,"Error parsing %s as start iwidth, final iwidth, ifactor\n",optarg);
	exit(-1);
      }
      break;
    case 'w':
      nscan = comma3floatscan(optarg,&wstart,&wend,&wfactor);
      if (nscan!=3){
	fprintf(stderr,"Error parsing %s as start width, final width, factor\n",optarg);
	exit(-1);
      }
      break;
    case 'r':
      nscan = sscanf(optarg,"%f",&rfi_factor);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as rfi factor\n",optarg);
	exit(-1);
      }
      break;
    case 'u':
      nscan = sscanf(optarg,"%f",&minlength);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as minlength in seconds\n",optarg);
	exit(-1);
      }
      break;
    case 't':
      tscrunch=1;
      nscan = sscanf(optarg,"%d",&tscrunch_factor);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as tscrunch factor\n",optarg);
	exit(-1);
      }
      break;
    case 'f':
      fscrunch = 1;
      nscan = sscanf(optarg,"%d",&fscrunch_factor);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as fscrunch_factor\n",optarg);
	exit(-1);
      }
      break;
    case 'g':
      hasplotdev = 1;
      plotdev = new char[strlen(optarg)+1];
      nscan = sscanf(optarg,"%s",plotdev);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as pgplot device\n",optarg);
	exit(-1);
      }
      break;
    case 's':
      hasstd = 1;
      nscan = sscanf(optarg,"%s",stdfile);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as standard filename\n",optarg);
	exit(-1);
      }
      break;
    case 'e':
      extension = 1;
      ext_chars = new char[strlen(optarg)+1];
      nscan = sscanf(optarg,"%s",ext_chars);
      ext_chars[strlen(optarg)]='\0';
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as extension\n",optarg);
	exit(-1);
      }
      break;
    case 'S':
      mkstd = 1;
      break;
    case 'X':
      Xplore = 1;
      break;
    case 'd':
      dedisperse = 1;
      break;
    case 'H':
      showHeader = 1;
      break;
    case 'Z':
      psrchiveZap = 1;
      break;
    case 'I':
      showIntHeader = 1;
      break;
    case 'E':
      nscan = sscanf(optarg,"%f",&edgeFrac);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as fraction of band to zap or nchans to zap\n",optarg);
	exit(-1);
      }
      edgeZap=1;
      break;
    case 'W':
      showWeightHeader = 1;
      break;
    case 'P':
      plot = 1;
      nscan = comma2intscan(optarg,&Splot_sub,&Splot_chan);
      if (nscan!=2){
	fprintf(stderr,"Error parsing %s as subint and channel to be plotted\n",optarg);
	exit(-1);
      }
      break;
    case 'Y':
      Yplot = 1;
      break;
    case 'G':
      Gplot = 1;
      break;
    case 'M':
      nscan=sscanf(optarg,"%f",&maxMeridianAngle);
      if (nscan!=1) {
	fprintf(stderr,"Error parsing %s and maximum meridian angle in deg\n",optarg);
	exit(-1);
      }
      maxMeridian = 1;
      break;
    case 'T':
      Tscrunch = 1;
      break;
    case 'U':
      getflux = 1;
      break;
    case 'o':
      gettoas = 1;
      break;
    case 'p':
      pscrunch = 1;
      break;
    case 'b':
      bscrunch = 1;
      nscan = sscanf(optarg,"%d",&bfactor);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as bscrunch factor\n",optarg);
	exit(-1);
      }
      break;
    case 'x':
      nscan = sscanf(optarg,"%f",&worst_SNR);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as worst SNR for residual\n",optarg);
	exit(-1);
      }
      break;
    case 'D':
      newdm = 1;
      nscan = sscanf(optarg,"%lf",&input_dm);
      if (nscan!=1){
	fprintf(stderr,"Error parsing %s as new DM\n",optarg);
	exit(-1);
      }
      break;
    case 'F':
      Fscrunch = 1;
      break;
    case 'K':
      ReWeight = 1;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
      help = 1;
      if (help) {
	printf(usage,argv[0]);
	printf("\n");
	printf(" %s is a program to manipulate observations of a pulsar\n",argv[0]);
	printf(" which mirrors the functionality of the psrchive routines:\n");
	printf(" pam - the archive manipulator\n");
	printf(" pav/psrplot - the archive plotter of profiles and grayscales\n");
	printf(" paz - the archive RFI zapper\n");
        printf(" pas - the standard profile creator\n");
	printf(" psrflux - the pulsar flux computer\n");
	printf(" pat - the toa calculator\n");
        printf(" pdmp - the search for optimum period and dm in an observation\n");
	printf("\n");
	printf(" %s assumes that you have an nbin * npol * nchan * nsubint array\n",
	       argv[0]);
	printf(" of data that attempted to be pulse phase coherent\n");
	printf(" and is built for speed not nanosecond accuracy\n");
	printf(" It uses the environment variable OMP_NUM_THREADS\n");
	printf(" to parallelise operations.\n\n Usage: %s options files\n\n",
	       argv[0]);
	printf(" -0 don't seek shifts when exploring residuals - much faster\n\n");
	printf(" -4 file1 rot1,zoom1 file2 rot2,zoom2 file3 rot3,zoom3 file4 rot4,zoom4\nPlots four file freq vs phase plot\n");
	printf(" -a print to screen amplitudes of 0,0,0th profile\n");
	printf(" -b bscrunch factor (power of 2 preferred)\n");
	printf(" -c rotate the profiles by 1/2 a turn - warning ruins TOAs!!!!\n" );
	printf(" -d dedisperse the archive to a weighted centre freq\n");
	printf(" -e ext at the end of manipulations create a new file with ext extension\n");
	printf(" -f fscrunch_factor (must divide nchan)\n");
	printf(" -g dev (pgplot graphics device to use eg 1/xs) \n");
	printf(" -h this help message\n");
	printf(" -i istart,iend,ifactor  - RFI detection plan - defaults 0.01,0.7,1.5\n");
	printf(" -k - accept a click after each plot instead of returning automatically\n");
	printf(" -o generate toas for each subint and frequency channel - requires a std profile\n");
	printf(" -p pscrunch to Stokes I\n");
	printf(" -s standard profile filename\n");
	printf(" -t tscrunch_factor (untested but runs - doesn't update az,zen,para)\n");
	printf(" -u min_subint_length(s) - eliminates subints with t< value \n");
	printf(" -v verbose mode\n");
	printf(" -w \"minwidth,maxwidth,f\" when searching for peak S/N \n");
	printf("   minbins and maxbins are the minimum and max # of bins to search \n");
	printf("   between with steps a factor f between them \n");
	printf("   if f=1, bins increment by 1 per trial \n");
	printf(" -x worst_snr when combined with -z deletes profiles with\n");
	printf("   a residual S/N is worse than worst_snr\n");
	printf(" -z tries to clean up RFI by comparing profile to std \n");
	printf("     \n");
	printf(" -A new_archive_name - incoherently adds archives together\n");
	printf(" -B compute the baseline from the std deviation of the baselines\n");
	printf(" -C \"Gain,Trec\" - provide the Gain and Receiver Temperature\n");
	printf(" -D new_DM - sets a new DM for the archive\n");
	printf(" -E fraction or -E nchans - delete the edges of the band\n");
	printf(" -F completely Fscrunch\n");
	printf(" -G grayscale of frequency vs phase\n");
	printf(" -H show the 80-byte header summary\n");
	printf(" -I show the subintegration header\n");
	printf(" -K reset weights to B(Mhz)*t(s)*npol\n");
	printf(" -L low S/N regime for RFI rejection\n");
	printf(" -M max_meridian_angle - trims subints great than max_meridian_angle from Observation\n");
	printf(" -N \"nx,ny\" for each archive plot grid of results\n");
	printf(" -O \"npturns,nfreqturns \" optimise for period and DM out to a max value\n");
	printf(" -P \"subint,channel\" plots profile of subint, channel");
	printf(" -Q psrbaselinefrac - default is 0.3\n");
	printf(" -R baseline - width of baseline when computing s/n's of the RFI\n");
	printf(" -S enter the standard profile creator\n");
	printf(" -T completely tscrunch\n");
	printf(" -U determine the flUx for every profile in the archive\n");
	printf(" -X explore your observation with the aid of a GUI - requires -z and -s\n");
	printf(" -W for every single profile show its weight and frequency\n");
	printf(" -Y plot time vs pulsar phase in a grayscale\n");
	printf(" -Z output a list of channels and subints for zapping in psh format\n");
      }
      exit(-1);
      break;
    case 'z':
      zap = 1;
      break;
    }
  }
  
  // open plot device in nx by ny cells unless Xploring
  if (!hasplotdev){
    plotdev = new char[2];
    strcpy(plotdev,"?");
  }
  if ((plot || multiplot || Splot || optimise || Yplot || Gplot) && !Xplore && !mkstd) {
    cpgbeg(0,plotdev,nxplot,nyplot);
  }
  if (Xplore || mkstd){
    cpgbeg(0,plotdev,1,1);
  }

  if (addarchs){
    grandcount=0;
  }

  for (int index=optind;index<argc;index++){
    filename=argv[index];

  if (verbose) fprintf(stderr,"Loading file %s...",filename);
  archive * a = new archive(filename);
  if (verbose) fprintf(stderr," loaded!\n");

  if (newpsr) a->update_source(newpsrname);
  if (ReWeight) a->ReWeight();
  if (edgeZap) a->edgeZap(edgeFrac);
  if (maxMeridian) a->maxMeridianAngle(maxMeridianAngle);

  archive * std = NULL;
  int maxbase=0;
  if (hasstd) {
    std = new archive(stdfile);
    maxbase = std->subints[0]->polnprofiles[0]->profiles[0]->maxbaseline();
    if (verbose) fprintf(stderr,"total bins in baseline of std profile 0,0 is %d\n",maxbase);
    if (maxbase==(int)std->nbin || maxbase<=1){
      fprintf(stderr,"Error - not enough baseline in std profile %s\n",std->filename);
      exit(-1);
    }
  }
  if (showHeader) a->showHeader();
  if (showIntHeader) a->showIntHeader();
  if (showWeightHeader) a->showWeightHeader();

  int minbins=wstart*a->nbin; // trial width minimum
  int maxbins=wend*a->nbin;  // unsure what this is?
  int iminbins=wstart*a->nbin; // trial width minimum
  int imaxbins=wend*a->nbin;  // unsure what this is?
  float basewidth = basew;

  // for every single profile determine the S/N of the best square pulse
  // the rms of the baseline and the baseline height
  if (newdm) {
    a->set_DM(input_dm);
  }
  // Zapping now in three parts.
  // First remove the std profile to create residuals and record S/N of profile
  // Then model the bandpass using the stddev of the baseline
  // Finally zap according to the recipe from -r factor and -X S/N_residual_max
  // and the bandpasses
  // Option to explore data with -X at the end
  // -z zap option by default assumes -r 1.1 -X 8
  
  if (zap && hasstd){
    // make a copy of the pscrunched archive to work on
    archive * acopy = new archive(a);
    if (verbose) fprintf(stderr,"pscrunching archive copy\n");
    acopy->pscrunch(2);
    if (verbose) fprintf(stderr,"dedispersing archive copy\n");
    double refFreq=acopy->dedisperse(); //ready for the zero-shift option
    if (verbose) fprintf(stderr,"dedispersing std with %lf =refFreq\n",refFreq);
    double stdFreq=std->dedisperse(acopy->dm,refFreq);
    if (verbose) fprintf(stderr,"done dedispersing std with %lf =refFreq\n",refFreq);
    //create some space to store snrs and bandpasses
    size_t ncells = acopy->nsubint*acopy->nchan;
    if (verbose) fprintf(stderr,"main loop ncells is %zu\n",ncells);
    float * snr_residual = new float[ncells];
    float * snr_negidual = new float[ncells];
    float * measured_bandpass = new float[ncells];
    float * snr_prof= new float[ncells];
    float * rmses= new float[ncells];
    float * sigma= new float[ncells];
    float * sigma2= new float[ncells];
    float * fitted_bandpass= new float[acopy->nchan];
    if (verbose) fprintf(stderr,"main loop starting\n");
    // this is the loop that takes all of the time I would think
    for (size_t i=0;i<(unsigned long int)acopy->nsubint;i++){
      if (verbose) fprintf(stderr,"finding residual of subint %zu of %zu\n",i,acopy->nsubint);
      #pragma omp parallel for
      for (size_t j=0;j<(unsigned long int)acopy->nchan;j++){
	rmses[i*acopy->nchan+j]=acopy->subints[i]->polnprofiles[j]->profiles[0]->rms();
	sigma[i*acopy->nchan+j]=acopy->subints[i]->polnprofiles[j]->profiles[0]->stddev();
	// this next line hurts if you search in shift space O(nbin^2)
	if (!lowsnr){
	snr_prof[i*acopy->nchan+j]=
	  acopy->subints[i]->polnprofiles[j]->profiles[0]->residual(
		  std->subints[0]->polnprofiles[0]->profiles[0],
		  seekflag,&sigma2[i*acopy->nchan+j]);
	}else{
	  snr_prof[i*acopy->nchan+j]=1.0;
	}
      }
    }
    //basewidth = float(maxbase)/float(acopy->nbin); // try a mild baseline width to compute S/N or residual
    // this loop is parallelizable but is fairly quick
    for (size_t i=0;i<acopy->nsubint;i++){
      if (verbose) fprintf(stderr,"Computing S/N of residual for subint %zu of %zu\n",i,acopy->nsubint);
      #pragma omp parallel for
      for (unsigned long int j=0;j<(unsigned long int)acopy->nchan;j++){
	int width, istart;
	float baseline;
	snr_residual[i*acopy->nchan+j] = 
	               acopy->subints[i]->polnprofiles[j]->profiles[0]->snr_vwauto(
						&width, &istart, 
						&measured_bandpass[i*acopy->nchan+j], 
						ifactor, ibasewidth, &baseline, 
							iminbins, imaxbins, irms);
	acopy->subints[i]->polnprofiles[j]->profiles[0]->scale(-1.0); // invert
	float rms2;
	int width2, istart2;
	float baseline2;
	snr_negidual[i*acopy->nchan+j] = 
                       acopy->subints[i]->polnprofiles[j]->profiles[0]->snr_vwauto(
						 &width2, &istart2, &rms2, 
						 ifactor, ibasewidth, &baseline2, 
						 iminbins, imaxbins, irms);
      }
    }
    if (verbose) fprintf(stderr,"determining the bandpass for %s\n",a->filename);
    // fixbandpass(measured_bandpass,acopy->nchan,acopy->nsubint,11,1.1,0.25,fitted_bandpass);
    // this is where things are deleted by setting the weight to zero
    fixbandpass(sigma2,acopy->nchan,acopy->nsubint,11,1.1,0.25,fitted_bandpass);
    int rficount=0;
    for (unsigned long int i=0;i<(unsigned long int)acopy->nsubint;i++)
      for (unsigned long int j=0;j<(unsigned long int)acopy->nchan;j++){
	size_t index = i*acopy->nchan + j;
	if (sigma2[index]>fitted_bandpass[j]*rfi_factor || 
	    snr_residual[index]>worst_SNR || snr_negidual[index]>worst_SNR){
	  a->subints[i]->polnprofiles[j]->set_weight(0.0);
	  rficount++;
	}
      }
    if (verbose) fprintf(stderr,"Deleted %d points of %zu = %f %\n",
			 rficount,acopy->nsubint*acopy->nchan,
			 100.0*rficount/float(ncells));
    if (zap && strncmp(a->telescope,"PKS",3)==0){
      fprintf(stderr,"This is the Parkes UWL zapping method\n");
      // Determine the median bandpass and the mad of the bandpass
      float PKSmed=give_median(int(ncells),sigma2);
      float PKSmad=mad(int(ncells),sigma2);
      // If the number of points beyond 2 * mad exceeds a threshold delete
      // the whole channel
      //fprintf(stderr,"median bandpass is %f, mad is %f\n",PKSmed,PKSmad);
      float * fminus2 = new float[acopy->nchan];
      float * fplus2 = new float[acopy->nchan];
      for (size_t i=0;i<acopy->nchan;i++){
	fminus2[i]=fstrideabove(sigma2,acopy->nchan,acopy->nsubint,i,PKSmed-2*PKSmad);
	fplus2[i]=fstrideabove(sigma2,acopy->nchan,acopy->nsubint,i,PKSmed+2*PKSmad);
	if (fminus2[i]<0.82 || fplus2[i] > 0.13) {
	  for (size_t j=0;j<a->nsubint;j++)
	    a->subints[j]->polnprofiles[i]->set_weight(0.0);
	}
      }
    }
    if (Xplore && hasstd) {
      float ** metadata = new float * [7];
      metadata[0]=snr_prof;
      metadata[1]=snr_residual;
      metadata[2]=snr_negidual;
      metadata[3]=measured_bandpass;
      metadata[4]=rmses;
      metadata[5]=sigma;
      metadata[6]=sigma2;
      char header[]="S/N_Profile S/N_Residual S/N_Negidual Bandpass rms std.dev std.dev2";
      a->pscrunch(2);
      a->explore(header,metadata,std->subints[0]->polnprofiles[0]->profiles[0]);
    }
  }
  if (zap && !hasstd){
    if (verbose) fprintf(stderr,"Zapping based on S/N and bandpass of off-pulse\n");
    archive * acopy = new archive(a);
    acopy->pscrunch(2);
    size_t ncells = acopy->nsubint*acopy->nchan;
    float sigma[ncells];
    float measured_bandpass[ncells];
    float snr_prof[ncells];
    float fitted_bandpass[acopy->nchan];
    // Compute the std deviation of the profile
    for (size_t i=0;i<(unsigned long int)acopy->nsubint;i++){
      #pragma omp parallel for
      for (size_t j=0;j<(unsigned long int)acopy->nchan;j++)
	sigma[i*acopy->nchan+j]=acopy->subints[i]->polnprofiles[j]->profiles[0]->stddev();
    }
    // Compute the S/N of every profile according to current recipe
    // this loop is parallelizable but is fairly quick
    for (size_t i=0;i<acopy->nsubint;i++){
      #pragma omp parallel for
      for (unsigned long int j=0;j<(unsigned long int)acopy->nchan;j++){
	int width, istart;
	float baseline;
	snr_prof[i*acopy->nchan+j] = 
	               acopy->subints[i]->polnprofiles[j]->profiles[0]->snr_vwauto(
						&width, &istart, 
						&measured_bandpass[i*acopy->nchan+j], 
						wfactor, psrbasewidth, &baseline, 
							minbins, maxbins, irms);
      }
    }
    // model the baseline from the std deviations
    if (verbose) fprintf(stderr,"determining the bandpass for %s\n",a->filename);
    // fixbandpass(measured_bandpass,acopy->nchan,acopy->nsubint,11,1.1,0.25,fitted_bandpass);
    // this is where things are deleted by setting the weight to zero
    fixbandpass(measured_bandpass,acopy->nchan,acopy->nsubint,11,1.1,0.25,fitted_bandpass);
    // zap all above -X value and -r value
    int rficount=0;
    for (unsigned long int i=0;i<(unsigned long int)acopy->nsubint;i++)
      for (unsigned long int j=0;j<(unsigned long int)acopy->nchan;j++){
	size_t index = i*acopy->nchan + j;
	if (measured_bandpass[index]>fitted_bandpass[j]*rfi_factor || 
	    snr_prof[index]>worst_SNR){
	  a->subints[i]->polnprofiles[j]->set_weight(0.0);
	  rficount++;
	}
      }
    if (verbose) fprintf(stderr,"Deleted %d points of %zu = %f\n",rficount,
			 acopy->nsubint*acopy->nchan,100.0*rficount/float(ncells));
    // Xplore result if desired
    if (Xplore){
      float ** metadata = new float * [3];
      metadata[0]=snr_prof;
      metadata[1]=measured_bandpass;
      metadata[2]=sigma;
      char header[]="S/N_Profile Bandpass std.dev";
      a->explore(header,metadata);
    } 
  }   
  if (psrchiveZap){
     printf("#!/usr/bin/env psrsh\n");
     for (size_t i=0;i<a->nsubint;i++)
       for (size_t j=0;j<a->nchan;j++){
         if (a->subints[i]->polnprofiles[j]->profiles[0]->weight==0.0) 
	   printf("zap such %d,%d\n",(int)i,(int)j);
       }
  }
  if (pscrunch) {
    if (verbose) fprintf(stderr,"pscrunching %s\n",filename);
    a->pscrunch(2);
  }
  if (bscrunch) {
    if (verbose) fprintf(stderr,"bscrunching by a factor %d\n",bfactor);
    a->bscrunch(bfactor);
  }
  if (Tscrunch) {
    if (verbose) fprintf(stderr,"completely tscrunching %s\n",filename);
    tscrunch_factor=a->nsubint;
    tscrunch=1;
  }
  if (tscrunch) {
    if (verbose) fprintf(stderr,"tscrunching %s by a factor %d\n",filename,
			 tscrunch_factor);
    a->tscrunch(tscrunch_factor);
  }
  if (dedisperse){
    if (verbose) fprintf(stderr,"dedispersing %s\n",filename);
    a->dedisperse();
  }
  if (Fscrunch) {
    if (verbose) fprintf(stderr,"fscrunching %s\n",filename);
    a->fadd(0);
  }
  if (fscrunch) {
    if (verbose) fprintf(stderr,"fscrunching %s by a factor %d\n",filename,fscrunch_factor);
    a->fadd(fscrunch_factor);
  }
  
  if (getflux && hasstd){
    float rms_baseline;
    int ncells = a->nsubint*a->nchan*a->npol;
    float fluxes[ncells];
    float flux_errors[ncells];
    float noisefactors[ncells];
    float snrs[ncells];
    // for each subint and frequency channel report the S/N one line at a time
    if (std->nbin!=a->nbin) std->bscrunch(int(std->nbin)/int(a->nbin));
    for (size_t isub=0;isub<a->nsubint;isub++)
      for (size_t ichan=0;ichan<a->nchan;ichan++){
	  profile * testprof = new profile(a->subints[isub]->polnprofiles[ichan]->profiles[0]);
	  if (!seekflag) testprof->dedisperse(
				    std->subints[0]->polnprofiles[0]->profiles[0]->freq,
				    a->dm,a->subints[0]->folding_period);
          // here snr is the sum of the s/n of all the non-baseline bits / sqrt(Non)
	  float snr=testprof->residual(std->subints[0]->polnprofiles[0]->profiles[0],seekflag,
                      &rms_baseline);
	  float Tsky = tsky(a->source,testprof->freq);
	  int onbins = std->nbin-std->subints[0]->polnprofiles[0]->profiles[0]->numzero();
	  int offbins = std->nbin-std->subints[0]->polnprofiles[0]->profiles[0]->numzero();
	  float BNpt = testprof->weight*1e6;
	  float TrecPlusTCMB;
	  if (strncmp(a->telescope,"MeerKAT",7)==0){
	    TrecPlusTCMB = MeerKATGain * marisaSEFD(testprof->freq,64);
	  } else
	    TrecPlusTCMB = Trec;
          
	  noisefactors[ichan+isub*a->nchan]=rms_baseline*sqrt(BNpt)/(TrecPlusTCMB+Tsky-TCMB);
	  if (BNpt!=0.0){
	    fluxes[ichan+isub*a->nchan] = snr * (TrecPlusTCMB+Tsky-TCMB)/Gain/sqrt(BNpt/testprof->nbin)*sqrt(onbins)/
	      testprof->nbin*1000; // mJy
	    float noise_Jy = (TrecPlusTCMB+Tsky-TCMB)/Gain/
	      sqrt(0.91*BNpt/testprof->nbin);
	    flux_errors[ichan+isub*a->nchan] = noise_Jy * sqrt(1.0/onbins + 1.0/offbins) * onbins/testprof->nbin*1000;
	    if (verbose) fprintf(stderr,"rms noise in channel %zu is %f mJy Tsky %f BNpt %f G %f\n",
				 ichan, noise_Jy*1000.0, Tsky, BNpt, Gain);
	    a->subints[isub]->polnprofiles[ichan]->profiles[0]->scale(noise_Jy/rms_baseline);
	    if (strncmp(a->telescope,"MeerKAT",7)==0) {
	      fluxes[ichan+isub*a->nchan]/=sqrt(0.91);
	    }
	  }
	  else
	    fluxes[ichan+isub*a->nchan] = -1.0;
	  snrs[ichan+isub*a->nchan]=snr;
      }
    float median_noisefactor=give_median(ncells,noisefactors);
    for (size_t isub=0;isub<a->nsubint;isub++)
      for (size_t ichan=0;ichan<a->nchan;ichan++)
	  printf("flux %s %zu %zu %f %f +/- %f mJy %lf MHz fiddlefactor %f\n",
		 a->filename,isub,ichan,snrs[ichan+isub*a->nchan],
		 fluxes[ichan+isub*a->nchan]*noisefactors[ichan+isub*a->nchan]/median_noisefactor,
		 flux_errors[ichan+isub*a->nchan]*noisefactors[ichan+isub*a->nchan]/median_noisefactor,
		 a->subints[isub]->polnprofiles[ichan]->profiles[0]->freq,
		 noisefactors[ichan+isub*a->nchan]/median_noisefactor);
  }
  if (getflux && !hasstd){
    int ncells = a->nsubint*a->nchan*a->npol;
    float fluxes[ncells];
    float flux_errors[ncells];
    float noisefactors[ncells];
    float snrs[ncells];
    // for each subint and frequency channel report the S/N one line at a time
    for (size_t isub=0;isub<a->nsubint;isub++)
      for (size_t ichan=0;ichan<a->nchan;ichan++){
	  profile * testprof = new profile(a->subints[isub]->polnprofiles[ichan]->profiles[0]);
	  int bestwidth,istart; 
	  float rms,baseline;
	  float irms=0.0;
	  float snr=testprof->snr_vwauto(&bestwidth, &istart, &rms, 
                 wfactor, psrbasewidth, &baseline, 
		       int(wstart*testprof->nbin), 
		       int(wend*testprof->nbin), irms);
	  float Tsky = tsky(a->source,testprof->freq);  // includes CMB
	  int onbins = bestwidth;
	  int offbins = testprof->nbin - onbins;
	  float BNpt = testprof->weight*1e6;
	  fprintf(stderr,"Trec %f Tsky %f onbins %d Gain %f BNpt %f\n",
		  Trec,Tsky,onbins,Gain,BNpt);
	  float TrecPlusTCMB;
	  if (strncmp(a->telescope,"MeerKAT",7)==0){
	    TrecPlusTCMB = MeerKATGain * marisaSEFD(testprof->freq,64);
	  }
	  noisefactors[ichan+isub*a->nchan]=rms*sqrt(BNpt)/(TrecPlusTCMB+Tsky-TCMB);
	  if (BNpt!=0.0){
	    fluxes[ichan+isub*a->nchan] = snr*(TrecPlusTCMB+Tsky-TCMB)/Gain/
	      sqrt(BNpt)*sqrt(float(onbins))/sqrt(testprof->nbin-onbins)*1000.0;
	    flux_errors[ichan+isub*a->nchan] = (TrecPlusTCMB+Tsky-TCMB)/Gain/sqrt(BNpt/testprof->nbin)*
	      sqrt(1.0/onbins+1.0/offbins)*onbins/testprof->nbin*1000.0;
	    if (strncmp(a->telescope,"MeerKAT",7)==0) {
	      fluxes[ichan+isub*a->nchan]/=sqrt(0.91);
	      flux_errors[ichan+isub*a->nchan]/=sqrt(0.91);
	    }
	  }
	  else
	    fluxes[ichan+isub*a->nchan] = -1.0;
	  snrs[ichan+isub*a->nchan]=snr;
      }
  float median_noisefactor=give_median(ncells,noisefactors);
    for (size_t isub=0;isub<a->nsubint;isub++)
      for (size_t ichan=0;ichan<a->nchan;ichan++){
	  printf("FluxFor %s %zu %zu %f %f +/- %f mJy %lf MHz Gain %f\n", 
		 a->filename,isub,ichan,snrs[ichan+isub*a->nchan],
		 fluxes[ichan+isub*a->nchan]*noisefactors[ichan+isub*a->nchan]/median_noisefactor,
		 flux_errors[ichan+isub*a->nchan],
                 a->subints[isub]->polnprofiles[ichan]->profiles[0]->freq, Gain);
      }
  }

  if (addarchs){
    a->showHeader();
    if (grandcount==0){
      // copy into archive sum;
      if (verbose) fprintf(stderr,"copying first archive into grand archive\n");
      grand = new archive(a);
      grandcount++;
    } else 
      {
	if (verbose)fprintf(stderr,"incoherently adding next archive %d %s\n",grandcount,a->filename);
	grand->incoadd(a);
	grandcount++;
      }
  }

  // when it is all scrunched etc unload
  if (extension){
    if (verbose) fprintf(stderr,"Unloading new archive %s%s\n",a->filename,ext_chars);
    if (minlength!=0.0) a->trim(minlength);
    a->unload(a->filename,ext_chars);
  }
  if (!plot4 && (Nplot) &&(plot || multiplot || Splot || Yplot || Gplot || mkstd || optimise)) {
    // erase every Nx * Ny plots
    // advance to the next subplot
    int xcoord = plotcounter % nxplot +1;
    int ycoord = plotcounter / nxplot +1;
    cpgpanl(xcoord,ycoord);
    cpgeras();
    //if (plotcounter%(nxplot*nyplot)==0) cpgeras();
    plotcounter++;
    if (plotcounter%(nxplot*nyplot)==0) {
      plotcounter=0;
      cpgask(1);
    }
  }

  if (plot4){
    cpgsvp(0.1+float(plot4-1),0.3+float(plot4-1),0.15,0.85);    
    plot4++;
  }

  if (centre) a->rotateleft(float(0.5));
  if (multiplot) a->multiplot();
  if (plot) a->Splot(Splot_sub,Splot_chan, click);
  if (Gplot) a->Gplot(0,click);
  if (Yplot) a->Yplot(0,click);
  if (mkstd) a->mkstd();
  if (ascii) a->subints[0]->polnprofiles[0]->profiles[0]->ascii();
  if (gettoas && hasstd){
    for (size_t isub=0;isub<a->nsubint;isub++)
      for (size_t ichan=0;ichan<a->nchan;ichan++){
	// get the shift
	if (a->subints[isub]->polnprofiles[ichan]->profiles[0]->weight!=0.0){
	float shift, shifterr;
	get_shift(a->nbin,
		  std->subints[0]->polnprofiles[0]->profiles[0]->amps,
		  a->subints[isub]->polnprofiles[ichan]->profiles[0]->amps,
		  &shift, &shifterr);
	if (verbose) fprintf(stderr,"SHIFT PHASE %f\n",shift);
	char telescope[100];
	strcpy(telescope,a->telescope);
	tolower(telescope);
	// report the shift to stdout
	if (shifterr!=0.0)printf("%s$%d$%d$ %lf %21.15Lf %f %s\n",a->filename,
				 (int)isub,(int)ichan,
	       a->subints[isub]->polnprofiles[ichan]->profiles[0]->freq,
		    double(shift)*a->subints[isub]->folding_period/(double(86400.0))+
				 a->subints[isub]->mjd+correction(a->subints[isub]->mjd,telescope)
				 ,shifterr*a->subints[isub]->folding_period*1e6,telescope);
	}
      }
  }
  if (optimise){
    float graininess = 0.5;
    profile * bestprf;
    float bestSNR,bestP,bestDM;
    int bestw, beststart;
    if (hasstd) basew = float(maxbase)/float(a->nbin); // try a mild baseline width to compute S/N or residual
    if (!hasstd) basew = psrbasewidth;
    a->pdmp(a->nbin, a->nsubint, a->nchan, 
	    timeturns, dmturns,
	    wstart, wend, wfactor, basew, graininess,
	    &bestw, &bestSNR,
	    &bestP, &bestDM, bestprf, &beststart);
    if (verbose) fprintf(stderr,"Successful return from pdmp bestw %d bestSNR %f bestDM %f beststart %d\n",
			 bestw,bestSNR,bestDM,beststart);
  }
  }
  if (plot || multiplot || Splot || Yplot || Gplot) cpgend();
  if (addarchs){
    grand->unload(grandfilename);
  }
  return(0);
}

