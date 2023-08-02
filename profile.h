void fscrunch_subint(float * subint, int nbin, int nchan, float * profile);
void fscrunch_archive(float * subints, int nbin, int nchan, int nsub,
		      float * profiles);

void tscrunch(float * profiles, int nbin, int nsub, float * profile);
void tscrunch(float * profiles, int nbin, int nsub,
	      float * profile, float drift);
void scrunchslope(float * profile, int nbin, float * subints, int nsub, float dph);

float bestslope(float * profiles, int nbin, int nsub,
		float maxslope, float dslope, float factor,
		float baselinefrac, float * bestsnr, int minbins, int maxbins, 
		float irms);

float get_ref_freq(float frch1, float df, int nchan);
float dmshift_turns(float period, float ref_freq, float freq, float dm);

void dedisperse_subint(int nbin, int nchan, float * subint, float frch1, float df,
		       float dm, float period, float * dedispersed);
void dedisperse_archive(int nbin, int nchan, int nsub,
			float * subints, float frch1, float df,
			float dm, float period, float * dedispersed);

// graphics
void prfplot(int nbin, float * prf, int istart, int iwidth, float baseline);

// optional
void zeroprf(float * prf, int nbin);
int make_prf(float * prf, int nbin, float amp, float width, float centre);
int make_arch(float * subints, int nbin, int nsub, float amp, float width,
	      float centre, float slope);
double psr_normal();
