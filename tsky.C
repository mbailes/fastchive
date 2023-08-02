#include <stdio.h>
#include <math.h>
#include <string.h>

#define TCMB 2.7
#define spindex -2.6

// Uses Dunc's old code and psrcat to determine
// Tsky at a particular frequency
// If reused for the same pulsar it redetermines
// pulsar position and tsky with ugly system calls
// requires psrcat and tsky to be installed on the system
// Minimal error checking at this point

float tsky(char * psrname, float freq){
  char chars[1000];
  char command[1000];
  static char lastpsrname[20];
  static float lastTemp=-1.0;
  static float lastfreq=-1.0;

  if (strlen(psrname)>19){
    fprintf(stderr,"tsky: Pulsar name incompatible with psrcat %s\n",psrname);
    return(-1);
  }

  if (freq==0.0){
    fprintf(stderr,"tsky: Don't try to measure tsky at zero frequency %f\n",freq);
    return(-1);
  }

  // if it is the same pulsar no need to run psrcat
  // for position
  if (strcmp(lastpsrname,psrname)==0){
    // if it is the same freq just return the last value
    if (lastfreq==freq) return(lastTemp);
    float Tsky=lastTemp-TCMB;
    // scale to new frequency
    float Tnew=pow(freq/lastfreq,spindex)*Tsky+TCMB;
    lastfreq=freq;
    lastTemp=Tnew;
    return(Tnew);
  } else {
    strcpy(lastpsrname,psrname);
  }

  sprintf(command,"psrcat -all -nohead -c \" gl gb \" %s",psrname);
  
  FILE * fptr = popen(command,"r");
  fgets(chars, 50, fptr);
  pclose(fptr);

  int i;
  float gl,gb;
  int nscanned = sscanf(chars,"%d %f %f",&i,&gl,&gb);
  //printf("gl = %f gb = %f\n",gl,gb);

  if (nscanned!=3) {
    fprintf(stderr,"Error uses psrcat to determine l and b for %s\n",
	    psrname);
    // try again with msp file from Renee
    sprintf(command,"psrcat -all -nohead -db_file /fred/oz002/users/mbailes/msp_ephems.db -c \" gl gb \" %s",psrname);
    fptr = popen(command,"r");
    fgets(chars, 50, fptr);
    pclose(fptr);
    nscanned = sscanf(chars,"%d %f %f",&i,&gl,&gb);
    if (nscanned !=3) return(-1);
  }

  sprintf(command,"tsky %f %f %f",gl,gb,freq);
  fptr = popen(command,"r");
  float temp;
  fgets(chars,50,fptr);
  pclose(fptr);
  nscanned=sscanf(chars,"%f",&temp);
  if (nscanned!=1){
    fprintf(stderr,"tsky error determining temp for l=%f b=%f freq=%f\n",gl,gb,freq);
    return(-1);
  }
  lastTemp=temp;
  lastfreq=freq;
  
  //printf("%f\n",temp);
  return(temp);
}

/*  unit test code eg a.out J0437-4715 843
int main(int argc, char ** argv){
  float freq;
  float temp;

  if (argc==3){
    int nscanned=sscanf(argv[2],"%f",&freq);
    if (nscanned==1) temp=tsky(argv[1], freq);
    printf("%f\n",temp);
  }
}
*/
