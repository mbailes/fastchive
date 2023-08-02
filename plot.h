float themin(int, float *);
float themax(int, float *);
float themin5(int, float *);
float themax5(int, float *);
int closest(int n, float x, float y, float * xx, float * yy, float dx, float dy);
float medianstride(float * data, size_t nchan, size_t nsubint, size_t channel);
float madstride(float * data, size_t nchan, size_t nsubint, size_t channel);
float fstrideabove(float * data, size_t nchan, size_t nsubint, size_t chan, float threshold);
float give_median(int npts, float * data);
void sort(int npts, float * data);
float mad(int npts, float * data);
int local_medians(int nsub, int nlocal, float * inputdata, float * median);
float give_fraction(int npts, float * newdata, float f);
void fixbandpass(float * data, size_t nchan, size_t nsubint, int nlocal, float tol,
		 float fraction, float * ideal_bandpass);

