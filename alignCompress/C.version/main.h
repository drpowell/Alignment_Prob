
#include <stdlib.h>
#include <math.h>

#define DOUBLE double
#define DO_COUNTS
 
#define loge2       M_LN2                    // ln(2);
#define log2e       M_LOG2E                  // log_2(e)

#ifndef   __USE_ISOC99
// C99 standard includes log2 and exp2 functions

#define INFINITY    1E20

static inline DOUBLE log2(DOUBLE a) {
  return log(a) * log2e;
}

static inline DOUBLE exp2(DOUBLE a) {
  return exp(a * loge2);
}
#endif


#ifdef __FAST_MATH__
#define exp2 __pow2		/* This is provided by math.h.  It is only valid for exp2(x) when -1<x<1 */
#endif

// Functions in markov.c
void fit_markovN(int alphaSize, int order, int len, unsigned char *str, DOUBLE *res);



// -----------------------------------------------------------------------------------
 
static inline DOUBLE logplus(DOUBLE a, DOUBLE b) {
  //  return -log(exp(-a) + exp(-b));    But this would lose accuracy
  if (isinf(a)) return b;
  if (isinf(b)) return a;
  if(a>b) {DOUBLE t=b;b=a;a=t;};   // make b>a
  if (b>a+50) return a;          // Approx if b >> a
  return a-log2(1+exp2(a-b));
}


// -----------------------------------------------------------------------------------
