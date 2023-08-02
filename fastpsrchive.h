
float archmax(int nbin, int nchan, int nsub, float * arch);
float archmin(int nbin, int nchan, int nsub, float * arch);
void pgarch(float * arch, int nbin, int nchan, int nsub);
float prfmin(float * pr, int nbin);
float prfmax(float * pr, int nbin);
float profdist(float f1, float f2, float turn);
void fscrunch_subint(float * subint, int nbin, int nchan, float * profile);
void fscrunch_archive(float * subints, int nbin, int nchan, int nsub,
		      float * profiles);
void tscrunch(float * profiles, int nbin, int nsub, float * profile);
void tscrunch(float * profiles, int nbin, int nsub,
	      float * profile, float drift);
void scrunchslope(float * profile, int nbin, float * subints, int nsub, float dph);
int load_arch(char * fname, int * nbin, int * npol, int * nchan, int * nsub,
	      float ** data);
int unload_arch(char * fname, int nbin, int npol, int nchan, int nsub,
		float * data);
double psr_normal();
float find_baseline(float * pr, int nbin, int bwidth, int * startbase);
void remove_constant(float * pr, int nbin, float b);
float find_rms(float * pr, int nbin, int nwidth, int istart);
float find_rms(float * pr, int nbin);
float rms_section(float * pr, int nbin, int istart, int width);
float snringate(float * pr, int nbin, int iwidth, int igate, float rms);
float snr_w(float * pr, int nbin, int nwidth, int * istart, float rms);
float snr_vws(float * pr, int nbin, int * bestwidth, int * istart,
	      float factor, float rms, int minbins, int maxbins);
float snr_vw(float * pr, int nbin, int * bestwidth, int * istart,
	     float factor, float rms, int maxbins);
float snr_vwauto(float * pr, int nbin, int * bestwidth, int * istart,
		 float * rms, float factor, float baselinefrac, float * baseline, 
		 int minbins, int maxbins, float irms);
float bestslope(float * profiles, int nbin, int nsub,
		float maxslope, float dslope, float factor,
		float baselinefrac, float * bestsnr, int minbins, int maxbins, 
		float irms);
void prfplot(float * prf, int nbin);
int make_prf(float * prf, int nbin, float amp, float width, float centre);
int make_arch(float * subints, int nbin, int nsub, float amp, float width,
	      float centre, float slope);
void zeroprf(float * prf, int nbin);
float get_ref_freq(float frch1, float df, int nchan);
float dmshift_turns(float period, float ref_freq, float freq, float dm);
void rotate(int nbin, float * dat, int offset, float * rotated);
void dedisperse_subint(int nbin, int nchan, float * subint, float frch1, float df,
		       float dm, float period, float * dedispersed);
void dedisperse_archive(int nbin, int nchan, int nsub,
			float * subints, float frch1, float df,
			float dm, float period, float * dedispersed);

