/*  m-align
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 *
    This file is part of m-align.
*/



#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "markov.h"

// fit_markovN - fit an arbitrary order markov model to a sequence.
//               return results in 'res'.  Note that res is really a two
//               dimensional array ie.  DOUBLE res[len][alphaSize]
//  This first fits the model, then applies it AND it does NOT charge for the model cost
//  thus it is cheating in the compression sense.
//  But for our purposes both the null model and the alignment model would have the same
//  parameter costs, so we can ignore them.
//  We don't use an adaptive model cause we want the start of the sequence treated the same
//  as the ends

static int alphaSize, order;
static int *charCounts=NULL, *countTotal=NULL;

static int past2int(int alphaSize, unsigned char *str, int len)
{
  int i, res = 0;
  for (i=0; i<len; i++)
    res = (res*alphaSize)+str[i];
  return res;
}

void markov_init(int _alphaSize, int _order)
{
  alphaSize = _alphaSize;
  order = _order;
  assert(order>=-1 && "Bad order="+order);

  markov_free();

  if (order>=0) {
    int i;
    int num_charCounts = pow(alphaSize, order+1);
    int num_countTotal = pow(alphaSize, order);
    charCounts = calloc(num_charCounts, sizeof(int));
    countTotal = calloc(num_countTotal, sizeof(int));
    for (i=0; i<num_charCounts; i++) charCounts[i]=1;
    for (i=0; i<num_countTotal; i++) countTotal[i]=alphaSize;
  }

}

void markov_free()
{
  if (charCounts) free(charCounts);
  if (countTotal) free(countTotal);
  charCounts=NULL;
  countTotal=NULL;
}

// markov_fit - increment markov counts from str
//              This may be called multiple time to fit multiple strings
void markov_fit(int len, unsigned char *str)
{
  int i;
  if (order>=0)
    for (i=0; i<len; i++) {
      if (i>=order) {
	int index = past2int(alphaSize, str+i-order, order);
	charCounts[ index*alphaSize+str[i] ]++;
	countTotal[ index ]++;
      }
    }
}

void markov_load(char *fname) 
{
  int i;
  FILE *f = fopen(fname, "r");
  char buf[80];
  if (!f) {
    fprintf(stderr, "Unable to read file '%s'\n", fname);
    exit(1);
  }

  {
    int o, alpha;
    int r;
    fgets(buf, 80, f);
    r = sscanf(buf, "ALPHASIZE=%d ORDER=%d\n", &alpha, &o);
    if (r != 2) {
      fprintf(stderr,"Bad header line in markov counts file: %s\n", buf);
      exit(1);
    }
    if (alpha != alphaSize) {
      fprintf(stderr,"Bad alphabet size in markov counts file.  Compiled with=%d.  Read=%d\n",
              alphaSize, alpha);
      exit(1);
    }
    if (o != order) {
      fprintf(stderr,"Bad markov order in markov counts file.  Read=%d Expected=%d\n", o, order);
      exit(1);
    }
  }

  while (fgets(buf, 80, f)) {
    char *p = strchr(buf, ':');
    int i,index;

    if (!p) {
      fprintf(stderr,"Bad line format in markov counts file: %s\n", buf);
      exit(1);
    }

    if (p-buf != order+1) {
      fprintf(stderr,"Markov count line is of different length to order: %s\n",buf);
      exit(1);
    }

    *p=0;
    for (i=0; buf[i]; i++) buf[i] = char2int(buf[i]);
    index = past2int(alphaSize, buf, order+1);

    charCounts[index] = atoi(p+1);
  }
  fclose(f);

  // Now setup countTotal[]
  for (i=0; i<pow(alphaSize, order); i++) {
    int j;
    int total = 0;
    for (j=0; j<alphaSize; j++)
      total += charCounts[i*alphaSize+j];
    countTotal[i] = total;
  }
}

void markov_save(FILE *f)
{
  int i;
  char buf[80];
  
  if (!f) {
    fprintf(stderr, "Bad FILE*'\n");
    exit(1);
  }

  fprintf(f, "ALPHASIZE=%d ORDER=%d\n", alphaSize, order);

  if (order<0) return;

  for (i=0; i < pow(alphaSize, order+1); i++) {
    int j, t = i;
    for (j=order; j>=0; j--) {
      buf[j] = int2char(t % alphaSize);
      t = t / alphaSize;
    }
    buf[order+1]=0;
    fprintf(f, "%s: %d\n", buf, charCounts[i]);
  }

}


void markov_predict(int len, unsigned char *str, DOUBLE *res)
{
  int i;

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



// markov_entropy - entropy of the markov model
// Note: Any high markov model can be considered as a first order
//       model with a large number of states
// Entropy = Sum_(over all contexts) ( P(context) * 
//                      Sum_(over all alpha chars) (P(a|context) * -log P(a|context) )
// To calculate this we need the steady-state state probabilities P(context).
// The transition probabilities can be used to form a number of linear equations:
//    P(context1)  = Sum_(over all contexts) ( P(context1| context) )
// With the equation Sum_(over all contexts) ( P(context) ) = 1
// These linear equations can be solved using Gaussian elimination.
DOUBLE markov_entropy()
{
  void solve_matrix(int n, DOUBLE **a_s, DOUBLE *b_s);

  int n = pow(alphaSize, order);
  DOUBLE *a_s[n];
  DOUBLE b_s[n];
  int i;

  for (i=0; i<n; i++)
    a_s[i] = calloc(n, sizeof(DOUBLE));

  // Setup the matrix of linear equations
  for (i=0; i<n; i++) {
    int j;
    for (j=0; j<alphaSize; j++) {
      int toState;
      if (order <= 0)
	toState = 0;
      else 
	toState = (i%(n/alphaSize))*alphaSize + j;
      a_s[toState][i] = 1.0 * charCounts[i*alphaSize+j]/countTotal[i];
    }
  }

  for (i=0; i<n; i++) {
    a_s[i][i] -= 1;
    b_s[i] = 0;
  }

  // Equations are under-constrained.  We can delete any 1 and replace it with sum(probs) = 1
  for (i=0; i<n; i++)
    a_s[0][i] = 1;
  b_s[0] = 1;

  solve_matrix(n, a_s, b_s);


  {
    // b_s now contains the equlibrium probabilities for each context
    DOUBLE entropy = 0;
    for (i=0; i<n; i++) {
      int j;
      DOUBLE h = 0;
      for (j=0; j<alphaSize; j++) {
	DOUBLE p = 1.0 * charCounts[i*alphaSize+j]/countTotal[i];
	h += p * -log2(p);
      }
      entropy += h * b_s[i];
    }

    return entropy;
  }
}


// --------------------------------------------------------------------------------

void print_matrix(int n, DOUBLE **a_s, DOUBLE *b_s) 
{
  int i;
  for(i=0 ; i<n; i++) {
    int j;
    printf("| ");
    for(j=0 ; j<n ; j++) {
      printf("%5.2f ", a_s[i][j]);
      //printf("%s", a_s[i][j] ? "x" : " ");
    }
    printf("| %s ", (i == (n/2) ? "=" : " "));
    printf(" | %5.2f |\n", b_s[i]);
  }
  printf("\n");
}

// Solve a set of linear equations of the form [As] * [xs] = [Bs]
// Solution in the b_s matrix
void solve_matrix(int n, DOUBLE **a_s, DOUBLE *b_s)
{
  int DEBUG = 0;
  int i;

  if (DEBUG) print_matrix(n, a_s, b_s);

  for(i=0; i<n-1; i++) {
    // Find the row with the largest value in column[i]  (this will be the pivot)
    int max = i;
    int i2;
    for (i2=i+1; i2<n; i2++) {
      if (fabs(a_s[i2][i] > fabs(a_s[max][i])))
	max = i2;
    }

    //printf("Swapping row[%d]=%f with max row[%d]=%f\n",i,a_s[i][i],max,a_s[max][i]);
    // Swap row $max with row $i  (and swap in the b vector)
    { 
      DOUBLE *t, tt;
      t = a_s[i];
      a_s[i] = a_s[max];
      a_s[max] = t;
      tt = b_s[i];
      b_s[i] = b_s[max];
      b_s[max] = tt;
    }

    // row a_s[i] is what we are subtracting from everything

    for(i2=i+1 ; i2 < n; i2++) {
      int j;
      // row a_s[i2] is the row we are about to subtract from
      DOUBLE r = a_s[i2][i] / a_s[i][i];
      if (r == 0) continue;

      for (j=i; j<n; j++) {
	// j is the column
	a_s[i2][j] -= a_s[i][j] * r;
      }
      b_s[i2] -= b_s[i] * r;
    }
  }

  if (DEBUG) print_matrix(n,a_s,b_s);

  for (i=n-1; i>=0; i--) {
    int j;
    DOUBLE sum = 0;
    for (j=i+1 ; j< n; j++) {
      sum += b_s[j] * a_s[i][j];
      a_s[i][j] = 0;
    }

    if (a_s[i][i] == 0) {
      if (b_s[i] == 0) {
	printf("Equations are underconstrained.\n");
	break;
      } else {
	printf("Equations are inconsistent.\n");
	break;
      }
    }

    b_s[i] = (b_s[i] - sum) / a_s[i][i];
    a_s[i][i] = 1;
  }
  if (DEBUG) print_matrix(n,a_s,b_s);
}

