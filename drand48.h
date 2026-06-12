/*
 * Declarations for drand48 routines on windows platform
 */


#ifndef _INC_DRAND48
#define _INC_DRAND48
#endif
 

double drand48(void);
double erand48(unsigned short xi[3]);
long lrand48(void);
long nrand48(unsigned short xi[3]);
long mrand48(void);
long jrand48(unsigned short xi[3]);
void srand48(long seedval);
unsigned short *seed48(unsigned short seed16v[3]);
void lcong48(unsigned short param[7]);
