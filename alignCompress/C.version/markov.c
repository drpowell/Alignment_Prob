#include <stdio.h>
#include <assert.h>

#include "main.h"

static int past2int(int alphaSize, unsigned char *str, int len)
{
  int i, res = 0;
  for (i=0; i<len; i++)
    res = (res*alphaSize)+str[i];
  return res;
}

// fit_markovN - fit an arbitrary order markov model to a sequence.
//               return results in 'res'.  Note that res is really a two
//               dimensional array ie.  DOUBLE res[len][alphaSize]
//  This first fits the model, then applies it AND it does NOT charge for the model cost
//  thus it is cheating in the compression sense.
//  But for our purposes both the null model and the alignment model would have the same
//  parameter costs, so we can ignore them.
//  We don't use an adaptive model cause we want the start of the sequence treated the same
//  as the ends
void fit_markovN(int alphaSize, int order, int len, unsigned char *str, DOUBLE *res)
{
  int i;
  int *charCounts=NULL, *countTotal=NULL;

  assert(order>=-1 && "Bad order="+order);

  if (order>=0) {
    int num_charCounts = pow(alphaSize, order+1);
    int num_countTotal = pow(alphaSize, order);
    charCounts = calloc(num_charCounts, sizeof(int));
    countTotal = calloc(num_countTotal, sizeof(int));
    for (i=0; i<num_charCounts; i++) charCounts[i]=1;
    for (i=0; i<num_countTotal; i++) countTotal[i]=alphaSize;
  }

  // Fit the model first
  if (order>=0)
    for (i=0; i<len; i++) {
      if (i>=order) {
	int index = past2int(alphaSize, str+i-order, order);
	charCounts[ index*alphaSize+str[i] ]++;
	countTotal[ index ]++;
      }
    }

  // Now cost per char.
  for (i=0; i<len; i++) {
    int j;
    int index=0;

    if (i>=order) 
      index = past2int(alphaSize, str+i-order, order);

    for (j=0; j<alphaSize; j++) {

      if (i<order || order<0) {
	res[i*alphaSize + j] = log2(alphaSize);
	continue;
      }


      res[i*alphaSize + j] = -log2(1.0*charCounts[index*alphaSize+j]/countTotal[index]);
    }

  }
}
