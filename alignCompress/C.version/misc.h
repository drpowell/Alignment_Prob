
#include <stdlib.h>
#include <math.h>

#define loge2  0.69314718055994530941   // log(2);
#define log2e  1.44269504088896340737   // 1.0/loge2;
 
inline double log2(double a) {
  //return Math.log(a) / loge2;
  return log(a) * log2e;
}
 
inline double exp2(double a) {
  //return Math.exp(a);
  return exp(a * loge2);
}
 
inline double logplus(double a, double b) {
  //  return -log(exp(-a) + exp(-b));    But this would lose accuracy
  if (isinf(a)) return b;
  if (isinf(b)) return a;
  if(a>b) {double t=b;b=a;a=t;};   // make b>a
  if (b>a+50) return a;          // Approx if b >> a
  return a-log2(1+exp2(a-b));
  //return a-Math.log(1+Math.exp(a-b));
}

