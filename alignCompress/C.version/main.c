// This is a re-implementation of the java version.
// It is supposed to be fast.  It is designed for the specific
// case of linear gaps, local alignment, summed over all alignments.



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "main.h"

#define DEBUG         0

#define ALPHA_SIZE    4
static char alphabet[ALPHA_SIZE] = "atgc";
//#define ALPHA_SIZE    24
//static char alphabet[ALPHA_SIZE] = "arndcqeghilkmfpstwyvbzx*";

char int2char(int i)
{
  if (i<0 || i>ALPHA_SIZE) {
    fprintf(stderr, "Bad param int2char(%d)\n", i);
    exit(-2);
  }
  return alphabet[i];
}

inline int char2int(char c)
{
  char *p = strchr(alphabet, tolower(c));
  if (!p) {
    fprintf(stderr, "Unknown char: %c\n", c);
    exit(-2);
  }
  return p - alphabet;
}

// -----------------------------------------------------------------------------------

static DOUBLE factorial(int N) {
  DOUBLE res = 1;
  int i;
  assert(N>=0);
  for (i=2; i<=N; i++) res *= i;
  return res;
}

// MMLparameter_cost - calculates the MML parameter cost for a Multinomial
  	// 'p' is the paramaters of the multinomial distribution.
	// 'N' is the number of things in this multinomial.
static DOUBLE MMLparameter_cost(DOUBLE *p, int numProbs, DOUBLE N) 
{
  DOUBLE h = factorial(numProbs-1); // h(theta) - prior probabilty density = (K-1)!
  DOUBLE F = 1.0/p[0];	// F will be the Fischer = N^(K-1)/(p1*p2*...*pk)
  DOUBLE cost;
  int i;
  for (i=1; i < numProbs-1; i++) {
    F *= N / p[i];
  }

  cost = 0.5 * log2(1 + F/(h*h*pow(12, numProbs-1)));
  
  cost += 0.5 * (numProbs-1) * log2e;
  
  return cost;
}


// -----------------------------------------------------------------------------------

// strict_DNA_seq - convert a DNA sequence into a sequence of only A,T,G or C.
//                  We do this by randomly choosing a character represented by the code
//  This is far from the best strategy, but probably sufficient if they are few and far between
void strict_DNA_seq(char *str, int len)
{
  int i;
  srandom(time(0));
  for (i=0; i<len; i++) 
    switch (tolower(str[i])) {
    case 'a': case 't':
    case 'g': case 'c':
      break;
    case 'u': str[i] = 't'; break;
    case 'r': str[i] = ((char[]){'g', 'a'})[random()%2]; break;
    case 'y': str[i] = ((char[]){'t', 'c'})[random()%2]; break;
    case 'k': str[i] = ((char[]){'g', 't'})[random()%2]; break;
    case 'm': str[i] = ((char[]){'a', 'c'})[random()%2]; break;
    case 's': str[i] = ((char[]){'g', 'c'})[random()%2]; break;
    case 'w': str[i] = ((char[]){'a', 't'})[random()%2]; break;
    case 'b': str[i] = ((char[]){'g', 't', 'c'})[random()%3]; break;
    case 'd': str[i] = ((char[]){'g', 'a', 't'})[random()%3]; break;
    case 'h': str[i] = ((char[]){'a', 'c', 't'})[random()%3]; break;
    case 'v': str[i] = ((char[]){'g', 'c', 'a'})[random()%3]; break;
    case 'n': str[i] = ((char[]){'a', 'g', 'c', 't'})[random()%4]; break;
    default:
      fprintf(stderr, "Unknown DNA character '%c'\n", str[i]);
      exit(-1);
    }
}

// -----------------------------------------------------------------------------------


// Ok these are the FSM costs    (taken from java version)
DOUBLE match_cost  = 0.5476268258420525;
DOUBLE change_cost = 1.6626638258420525;

DOUBLE diag_fromD  = 0.16198884271367536;
DOUBLE start_fromD = 4.235036668555728;

DOUBLE diag_fromI  = 0.8317598880466708;
DOUBLE start_fromI = 4.904807713888723;
DOUBLE cont_fromI  = 1.3048077138887235;

// -----------------------------------------------------------------------------------
enum cIndex {match_i, change_i, diag_fromD_i, start_fromD_i, 
	     diag_fromI_i, start_fromI_i, cont_fromI_i, num_counts};

struct val_counts {
  DOUBLE val;
#ifdef DO_COUNTS
  DOUBLE counts[num_counts];
#endif
};

struct cell {
  struct val_counts d,h,v;
};

// -----------------------------------------------------------------------------------
// counts_to_params - convert counts into parameters.  
// Note: This modifies the global cost variabls!
//   Note the 0.5* in the start_fromD computation.  This is cause there are 2 start_fromD 
//   arcs out of each state, but we only have one count for them combined.
void counts_to_params(DOUBLE c[])
{
  {
    DOUBLE sum = c[match_i] + c[change_i];
    match_cost  = -log2(c[match_i]/sum);
    change_cost = -log2(c[change_i]/sum);
  }

  {
    DOUBLE sum1 = c[diag_fromD_i] + c[start_fromD_i];
    DOUBLE sum2 = c[diag_fromI_i] + c[start_fromI_i] + c[cont_fromI_i];
    diag_fromD   = -log2(c[diag_fromD_i]/sum1);
    start_fromD  = -log2(0.5 * c[start_fromD_i]/sum1);
    diag_fromI   = -log2(c[diag_fromI_i]/sum2);
    start_fromI  = -log2(c[start_fromI_i]/sum2);
    cont_fromI   = -log2(c[cont_fromI_i]/sum2);
  }
}

// -----------------------------------------------------------------------------------

void init_counts(DOUBLE c[])
{
  int i;
  for (i=0; i<num_counts; i++)
    c[i] = 0.5;
}

// combine does the following:
//    v1 = val1->val  c1 = val1->counts
//       Let p1 = 2**-v1   and p2 = 2**-v2
//     c1[i] = (p1 * c1[i] + p2 * c2[i]) * 1/(p1+p2)
// Then:
//   val1->val = logplus(val1->val, val2->val)
void combine(struct val_counts *val1, struct val_counts *val2)
{
  DOUBLE v1 = val1->val;
  DOUBLE v2 = val2->val;
#ifdef DO_COUNTS
  DOUBLE *c1 = val1->counts;
  DOUBLE *c2 = val2->counts;
#endif

  if (isinf(v2)) return;
  if (isinf(v1)) {
    val1->val = val2->val;
#ifdef DO_COUNTS
    memcpy(c1, c2, num_counts * sizeof(DOUBLE));
#endif
    return;
  }

#ifdef DO_COUNTS
  if (v1<v2) {
    DOUBLE scale = exp2(v1-v2);
    DOUBLE norm = 1.0/(1 + scale);
    int i;
    for (i=0; i<num_counts; i++)
      c1[i] = norm * (c1[i] + c2[i] * scale);
  } else {
    DOUBLE scale = exp2(v2-v1);
    DOUBLE norm = 1.0/(1 + scale);
    int i;
    for (i=0; i<num_counts; i++)
      c1[i] = norm * (scale * c1[i] + c2[i]);
  }
#endif

  val1->val = logplus(val1->val, val2->val);
}

// combine_inc - calls combine and does:
//           param = *val2
//           param.val += incVal;    param.counts[index]++
//         combine( val1, param)
void combine_inc(struct val_counts *val1, struct val_counts *val2, DOUBLE incVal, int index)
{
  DOUBLE v = val2->val;
  val2->val += incVal;
#ifdef DO_COUNTS
  if (index>=0) val2->counts[index]++;
#endif

  combine(val1, val2);
  
#ifdef DO_COUNTS
  if (index>=0) val2->counts[index]--;
#endif
  val2->val = v;
}
// -----------------------------------------------------------------------------------


DOUBLE calcModel(struct val_counts *v)
{
#ifndef DO_COUNTS
  return 0;
#else
  DOUBLE probs[3];
  DOUBLE dataLen;
  DOUBLE len = 0;

  // First match/change params
  probs[0] =  exp2(-match_cost);
  probs[1] = exp2(-change_cost);
  dataLen = v->counts[diag_fromD_i] + v->counts[diag_fromI_i];
  len += MMLparameter_cost(probs, 2, dataLen);

  // Now params from diag state
  probs[0] = exp2(-diag_fromD);
  probs[1] = 2 * exp2(-start_fromD);
  dataLen = v->counts[diag_fromD_i] + v->counts[start_fromD_i];
  len += MMLparameter_cost(probs, 2, dataLen);

  // Finally params from ins/del state
  probs[0] = exp2(-diag_fromI);
  probs[1] = exp2(-start_fromI);
  probs[2] = exp2(-cont_fromI);
  dataLen = v->counts[diag_fromI_i] + v->counts[start_fromI_i] + v->counts[start_fromI_i];
  len += MMLparameter_cost(probs, 3, dataLen);

  return len;
#endif
}

// -----------------------------------------------------------------------------------

DOUBLE doAlign(unsigned char *seqA, unsigned char *seqB, int lenA, int lenB,
	       DOUBLE seqA_enc[lenA][ALPHA_SIZE], DOUBLE seqA_cum[lenA+1],
	       DOUBLE seqB_enc[lenB][ALPHA_SIZE], DOUBLE seqB_cum[lenB+1],
	       DOUBLE *final_counts /* If the caller wants the final counts */
	       )
{
  struct cell D[2][lenB+1];
  int i,j;

  struct val_counts empty, final;

  DOUBLE mdlCost;
  
#ifdef DO_COUNTS
  init_counts(empty.counts);
  init_counts(final.counts);
#endif
  final.val = INFINITY;

  // Clear first row
  for (j=0; j<lenB+1; j++) {
    D[0][j] = ((struct cell) {.d = {.val = INFINITY}, 
			      .v = {.val = INFINITY},
			      .h = {.val = INFINITY}});
  }

  for (i=0; i<lenA+1; i++) {

    // Clear next row
    for (j=0; j<lenB+1; j++) {
      D[(i+1)%2][j] = ((struct cell) {.d = {.val = INFINITY}, 
				      .v = {.val = INFINITY},
				      .h = {.val = INFINITY}});
    }
    
    for (j=0; j<lenB+1; j++) {
      
      // Add the jump _into_ D[i][j]
      empty.val = seqA_cum[i]+seqB_cum[j];
      combine(&D[i%2][j].d, &empty);
      

      // D[i][j] is now fully computed. Do it's jump _out_ into final cell
      //    For simplicity, enforce jump to be out of diag state.  I think this is reasonable
      combine_inc(&final, &D[i%2][j].d, 
		  (seqA_cum[lenA]-seqA_cum[i]) + (seqB_cum[lenB]-seqB_cum[j]), -1);
      combine_inc(&final, &D[i%2][j].v, 
		  (seqA_cum[lenA]-seqA_cum[i]) + (seqB_cum[lenB]-seqB_cum[j]), -1);
      combine_inc(&final, &D[i%2][j].h, 
		  (seqA_cum[lenA]-seqA_cum[i]) + (seqB_cum[lenB]-seqB_cum[j]), -1);
      
      if (DEBUG>=3)
	printf("D[%d][%d] d=%f h=%f v=%f     final=%f\n", i,j, 
	       D[i%2][j].d.val, D[i%2][j].h.val, D[i%2][j].v.val, final.val);


      // And finally outputs of D[i][j] into the 3 cells that depend on it

      // Horizontal step
      if (j<lenB) {
	DOUBLE char_cost = seqB_enc[j][seqB[j]];
	combine_inc(&D[i%2][j+1].h, &D[i%2][j].d, start_fromD + char_cost, start_fromD_i);
	combine_inc(&D[i%2][j+1].h, &D[i%2][j].v, start_fromI + char_cost, start_fromI_i);
	combine_inc(&D[i%2][j+1].h, &D[i%2][j].h, cont_fromI  + char_cost, cont_fromI_i);
      }

      // Vertical step
      if (i<lenA) {
	DOUBLE char_cost = seqA_enc[i][seqA[i]];
	combine_inc(&D[(i+1)%2][j].v, &D[i%2][j].d, start_fromD + char_cost, start_fromD_i);
	combine_inc(&D[(i+1)%2][j].v, &D[i%2][j].v, cont_fromI  + char_cost, cont_fromI_i);
	combine_inc(&D[(i+1)%2][j].v, &D[i%2][j].h, start_fromI + char_cost, start_fromI_i);
      }

      // Diagonal step
      if (i<lenA  &&  j<lenB) {
	DOUBLE diagCost;
	int incIndex;
	if (seqA[i]==seqB[j]) {
	  // Do: P(match) * ( P(char a) + P(char b) ) / 2
	  diagCost = match_cost + logplus(seqA_enc[i][ seqA[i] ], seqB_enc[j][ seqB[j] ]) + 1;
	  incIndex = match_i;
	} else {
	  // Do: P(change) * P(char a) * P(char b) * 0.5 * (1/(1-P(char b)) + 1/(1-P(char a)))
	  DOUBLE aN = exp2(-seqA_enc[i][ seqB[j] ]);
	  DOUBLE bN = exp2(-seqB_enc[j][ seqA[i] ]);
	  DOUBLE norm = -log2( 1/(1-aN) + 1/(1-bN) );
	  diagCost = change_cost + seqA_enc[i][ seqA[i] ] + seqB_enc[j][ seqB[j] ] + 1 + norm;
	  incIndex = change_i;
	}

#ifdef DO_COUNTS
	D[i%2][j].d.counts[incIndex]++;
	D[i%2][j].v.counts[incIndex]++;
	D[i%2][j].h.counts[incIndex]++;
#endif	
	combine_inc(&D[(i+1)%2][j+1].d, &D[i%2][j].d, diag_fromD + diagCost, diag_fromD_i);
	combine_inc(&D[(i+1)%2][j+1].d, &D[i%2][j].v, diag_fromI + diagCost, diag_fromI_i);
	combine_inc(&D[(i+1)%2][j+1].d, &D[i%2][j].h, diag_fromI + diagCost, diag_fromI_i);
	
#ifdef DO_COUNTS
	D[i%2][j].d.counts[incIndex]--;
	D[i%2][j].v.counts[incIndex]--;
	D[i%2][j].h.counts[incIndex]--;
#endif
      }

    }
  }

  {
    // Add a cost for the length of the alignment.
    // Assume we know the length of the sequences. Need to encode the start and end
    // of the alignment.  Assume uniform over all positions.  Have 4 cut-points to
    // encode, but only need to encode 3 (kind of).
    // So encode 2 from the shorter sequence, and 1 from the longer.
    DOUBLE l1 = (lenA < lenB ? lenA : lenB);
    DOUBLE l2 = (lenA > lenB ? lenA : lenB);
    final.val += log2(l1) * 2 - (l1>1 ? 1 : 0);
    final.val += log2(l2);
  }

  mdlCost = calcModel(&final);
  printf("Aligment Encoding = %f bits    (model=%f data=%f)\n", mdlCost+final.val, mdlCost, final.val);

#ifdef DO_COUNTS
  if (final_counts) {		/* Does our caller want the final counts? */
    int i;
    for (i=0; i<num_counts; i++)
      final_counts[i] = final.counts[i];
  }

  if (DEBUG>=1) {		// Print the final counts
    int i;
    printf("Counts:\n");
    for (i=0; i<num_counts; i++)
      printf("  %d:%.3f", i, final.counts[i]);
    printf("\n");
  }
#endif

  return (seqA_cum[lenA]+seqB_cum[lenB]) - (mdlCost + final.val);
}


// -----------------------------------------------------------------------------------
DOUBLE alignDriver(int maxIterations, char *seqA, char *seqB)
{
  DOUBLE lastDiff;
  int i;
  int lenA = strlen(seqA);
  int lenB = strlen(seqB);
  unsigned char seqA_i[lenA];
  unsigned char seqB_i[lenB];


  DOUBLE seqA_enc[lenA][ALPHA_SIZE], seqB_enc[lenB][ALPHA_SIZE];
  DOUBLE seqA_cum[lenA+1], seqB_cum[lenB+1];

  // Convert DNA sequence to only an A G C or T
  strict_DNA_seq(seqA, lenA);
  strict_DNA_seq(seqB, lenB);

  if (DEBUG>=2) printf("Seq1:\n%s\nSeq2:\n%s\n", seqA, seqB);

  // First convert strings to numbers representing the characters
  for (i=0; i<lenA; i++) seqA_i[i] = char2int(seqA[i]);
  for (i=0; i<lenB; i++) seqB_i[i] = char2int(seqB[i]);

  fit_markovN(ALPHA_SIZE, 1, lenA, seqA_i, (DOUBLE*)seqA_enc);
  fit_markovN(ALPHA_SIZE, 1, lenB, seqB_i, (DOUBLE*)seqB_enc);
  

  // Calculate the cumulative encodings
  seqA_cum[0] = seqB_cum[0] = 0;
  for (i=0; i<lenA; i++) seqA_cum[i+1] = seqA_enc[i][seqA_i[i]] + seqA_cum[i];
  for (i=0; i<lenB; i++) seqB_cum[i+1] = seqB_enc[i][seqB_i[i]] + seqB_cum[i];


  // Print encoding for each character in input sequences.
  if (DEBUG>=2) {
    int j;
    for (j=0; j<ALPHA_SIZE; j++) {
      printf("%c ",int2char(j));
      for (i=0; i<lenA; i++) printf("%.2f ", seqA_enc[i][j]);
      printf("\n");
    }
    for (i=0; i<lenA+1; i++) printf("%.2f ", seqA_cum[i]);
    printf("\n\n");

    for (j=0; j<ALPHA_SIZE; j++) {
      printf("%c ",int2char(j));
      for (i=0; i<lenB; i++)
	printf("%.2f ", seqB_enc[i][j]);
      printf("\n");
    }
    for (i=0; i<lenB+1; i++) printf("%.2f ", seqB_cum[i]);
    printf("\n\n");
  }

  printf("Null mdl=%f bits (str1=%f bits  str2=%f bits)\n", 
	 seqA_cum[lenA]+seqB_cum[lenB], seqA_cum[lenA], seqB_cum[lenB]);

  lastDiff = -INFINITY;

  while (maxIterations--) {
    DOUBLE counts[num_counts];
    DOUBLE diff = doAlign(seqA_i, seqB_i, lenA, lenB,
			  seqA_enc, seqA_cum,
			  seqB_enc, seqB_cum,
			  counts);

    if (diff >= lastDiff && diff-lastDiff < 0.1) {
      lastDiff = diff;
      break;
    }

    counts_to_params(counts);

    lastDiff = diff;
  }

  return lastDiff;
}

// -----------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  DOUBLE r;
  char *seqA, *seqB;

  if (argc == 1) {
    seqA = malloc(10000);
    seqB = malloc(10000);
    fgets(seqA, 10000, stdin);
    fgets(seqB, 10000, stdin);    
    if (seqA[strlen(seqA)-1] != '\n') { fprintf(stderr, "seqA too long.\n"); exit(-1); }
    if (seqB[strlen(seqB)-1] != '\n') { fprintf(stderr, "seqB too long.\n"); exit(-1); }
    seqA[strlen(seqA)-1]=0;
    seqB[strlen(seqB)-1]=0;
  } else if (argc ==3) {
    seqA = argv[1];
    seqB = argv[2];
  } else {
    fprintf(stderr, "Usage: %s <seqA> <seqB>\n", argv[0]);
    exit(-1);
  }

  r = alignDriver(8, seqA, seqB);
  printf("log odds ratio = %f\n", r);

  return 0;
}

