/*  m-align
 *  Copyright (c) David Powell <david@drp.id.au>
 *
 * 
 *
    This file is part of m-align.




  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "common.h"
#include "markov.h"


#define DEBUG         0

int markovOrder = 0;
char *markovFile = NULL;
char *probsFileA = NULL;
char *probsFileB = NULL;
char *markovSaveFile = NULL;


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
  int i;
  DOUBLE cost;
  DOUBLE h = factorial(numProbs-1); // h(theta) - prior probabilty density = (K-1)!
  DOUBLE F = 1.0/p[0];	// F will be the Fischer = N^(K-1)/(p1*p2*...*pk)
  for (i=1; i < numProbs-1; i++) {
    F *= N / p[i];
  }

  cost = 0.5 * log2(1 + F/(h*h*pow(12, numProbs-1)));
  
  cost += 0.5 * (numProbs-1) * log2e;
  
  return cost;
}


// -----------------------------------------------------------------------------------


// Ok these are the FSM costs    (taken from java version)
//  They correspond to Smith-Waterman scores of: 5 for a match, -4 for a mismatch, 
//                 -16 and -4 for the initial and subsequence characters in a gap.
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
	       DOUBLE seqA_enc[][ALPHA_SIZE], DOUBLE seqA_cum[],
	       DOUBLE seqB_enc[][ALPHA_SIZE], DOUBLE seqB_cum[],
	       DOUBLE *final_counts /* If the caller wants the final counts */
	       )
{
  DOUBLE mdlCost;
  struct cell D[2][lenB+1];
  int i,j;

  struct val_counts empty, final;
  
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
  printf("Alignment Encoding = %f bits    (model=%f data=%f)\n", 
	 mdlCost+final.val, mdlCost, final.val);

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
void load_probs_file(char *fname, int len, void *r)
{
  FILE *f = fopen(fname, "r");
  char buf[1024];
  DOUBLE (*res)[len][ALPHA_SIZE];
  int linenum;
  int blank_lines;

  res = r;

  if (!f) {
    fprintf(stderr, "Unable to open file '%s'. Exiting...\n", fname);
    exit(1);
  }

  linenum = 0;
  blank_lines = 0;
  while (!feof(f)) {
    int l;
    char *p;
    int i;

    if (!fgets(buf, 1024, f))
      break;

    if (buf[0] == '#') {
      blank_lines++;
      continue;			// Comment line.  Skip it
    }

    l = strlen(buf);

    i = 0;
    for (p = strtok(buf, " \t"); p; p = strtok(NULL, " \t")) {
      int j;
      int blank_field = 1;
      for (j=0; p[j] && blank_field; j++)
	if (!isspace(p[j]))
	    blank_field = 0;
	    
      if (blank_field)
	continue;

      if (i>=ALPHA_SIZE) {
	fprintf(stderr, "Probs file '%s', line %d: too many probabilities.  Got %d, expect %d. Exiting...\n", 
		fname, linenum + blank_lines + 1, i+1, ALPHA_SIZE);
	exit(1);
      }

      if (linenum>=len) {
	fprintf(stderr, "Probs file '%s': too many lines.  Exiting...\n", fname);
	exit(1);
      }

      (*res)[linenum][i] = atof(p);
      i++;
    }

    if (i==0) {
      blank_lines++;
      continue;			// Blank line.  Skip it
    }

    if (i != ALPHA_SIZE) {
      fprintf(stderr, "Probs file '%s', line %d: too few probabilities.  Exiting...\n", 
	      fname, linenum + blank_lines + 1);
      exit(1);
    }


    { // Check probs sum to 1.0  (or there abouts).  And convert then to encoding lengths
      double sum=0;
      for (i=0; i<ALPHA_SIZE; i++) {
	sum += (*res)[linenum][i];
	(*res)[linenum][i] = -log2((*res)[linenum][i]);
      }
      if (fabs(1.0 - sum) > 0.001) {
	fprintf(stderr, "Probs file '%s', line %d: Probabilities don't sum to 1.  Sum=%f. Exiting...\n", 
		fname, linenum + blank_lines + 1, sum);
	exit(1);
      }
    }

    linenum++;
  }

  if (linenum != len) {
    fprintf(stderr, "Probs file '%s': too few lines. Got %d, expect %d. Exiting...\n",
	    fname, linenum, len);
    exit(1);
  }
}

// -----------------------------------------------------------------------------------
DOUBLE alignDriver(int maxIterations, char *seqA, char *seqB)
{
  int i;
  int lenA = strlen(seqA);
  int lenB = strlen(seqB);
  unsigned char seqA_i[lenA];
  unsigned char seqB_i[lenB];
  DOUBLE lastDiff;


  DOUBLE seqA_enc[lenA][ALPHA_SIZE], seqB_enc[lenB][ALPHA_SIZE];
  DOUBLE seqA_cum[lenA+1], seqB_cum[lenB+1];

  // Convert DNA sequence to only an A G C or T
  strict_DNA_seq(seqA, lenA);
  strict_DNA_seq(seqB, lenB);

  if (DEBUG>=2) printf("Seq1:\n%s\nSeq2:\n%s\n", seqA, seqB);

  // First convert strings to numbers representing the characters
  for (i=0; i<lenA; i++) seqA_i[i] = char2int(seqA[i]);
  for (i=0; i<lenB; i++) seqB_i[i] = char2int(seqB[i]);



  if (!probsFileA || !probsFileB) {
    markov_init(ALPHA_SIZE, markovOrder);
    if (markovFile)
      markov_load(markovFile);
  }

  if (!probsFileA && !markovFile)
    markov_fit(lenA, seqA_i);	/* One model for both sequences */

  if (!probsFileB && !markovFile)
    markov_fit(lenB, seqB_i);

  if (DEBUG>=1 && (!probsFileA || !probsFileB))
    printf("Markov model entropy = %f bits/char\n", markov_entropy());
  


  if (probsFileA)
    load_probs_file(probsFileA, lenA, seqA_enc);
  else
    markov_predict(lenA, seqA_i, (DOUBLE*)seqA_enc);

  if (probsFileB)
    load_probs_file(probsFileB, lenB, seqB_enc);
  else
    markov_predict(lenB, seqB_i, (DOUBLE*)seqB_enc);

  if (markovSaveFile) {
    FILE *f = fopen(markovSaveFile, "w");
    if (!f) {
      fprintf(stderr, "Unable to open file '%s' for writing.\n", markovSaveFile);
    } else {
      fprintf(stderr, "Saving Markov Model parameters to file '%s'\n", markovSaveFile);
      markov_save(f);
      fclose(f);
    }
  }

  markov_free();




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


void cite()
{
  fprintf(stdout, "\n");
  fprintf(stdout, "  Copyright (c) 2002 David Powell\n");
  fprintf(stdout, "  m-align aligns two sequences using the M-alignment\n");
  fprintf(stdout, "  algorithm.  (Version %s)\n", VERSION);
  fprintf(stdout, "  m-align comes with ABSOLUTELY NO WARRANTY; and is provided\n");
  fprintf(stdout, "  under the GNU Public License v2, for details see file COPYRIGHT\n\n");

  fprintf(stdout, "Please cite:\n");
  fprintf(stdout, "  L. Allison, D. R. Powell and T. I. Dix\n");
  fprintf(stdout, "  \"Compression and Approximate Matching\"\n");
  fprintf(stdout, "  The Computer Journal, 1999 42:1 pp1-10\n");
  fprintf(stdout, "\n");
}

void usage(char *exeName)
{
  fprintf(stderr, "Usage: %s [-mfiabh] <seqA> <seqB>\n", exeName);
  fprintf(stderr, "  Where <seqA> and <seqB> are files in FASTA format.\n");
  fprintf(stderr, "       -m <n>     Use a n'th order Markov Model (default=0)\n");
  fprintf(stderr, "       -f <file>  Read the Markov Model parameters from <file>\n");
  fprintf(stderr, "       -s <file>  Save the Markov Model parameters to <file> after fitting the sequences\n");
  fprintf(stderr, "       -i <n>     Set max iterations to <n>.  Use -1 for unlimited\n");
  fprintf(stderr, "       -a <file>  Read character prediction probabilites for seqA from <file>\n");
  fprintf(stderr, "       -b <file>  Read character prediction probabilites for seqB from <file>\n");
  fprintf(stderr, "       -h         This help\n");
  exit(-1);
}

int main(int argc, char **argv)
{
  char *exeName = argv[0];
  DOUBLE r;
  char *seqA, *seqB;
  int maxIterations = -1;

  cite();

  while (1) {
    int c = getopt(argc, argv, "m:f:s:i:a:b:h");
    if (c==-1)
      break;
    switch (c) {
    case 'm':
      markovOrder = atoi(optarg);
      break;
    case 'f':
      markovFile = optarg;
      break;
    case 's':
      markovSaveFile = optarg;
      break;
    case 'i':
      maxIterations = atoi(optarg);
      break;
    case 'a':
      probsFileA = optarg;
      break;
    case 'b':
      probsFileB = optarg;
      break;
    case 'h':
    default:
      usage(exeName);
    }
  }
	
  argc -= optind-1;
  argv += optind-1;

  if (argc != 3) {
    usage(exeName);
  } else {
    seqA = read_fasta(argv[1]);
    seqB = read_fasta(argv[2]);
  } 

  fprintf(stderr, "Markov Order = %d\n", markovOrder);
  if (markovFile)
    fprintf(stderr, "Using markov counts file '%s'\n", markovFile);

  fprintf(stderr, "Max iterations = %d\n", maxIterations);

  r = alignDriver(maxIterations, seqA, seqB);
  printf("\n\nlog odds ratio = %f bits (%s)\n", r,
      (r<0 ? "unrelated" : "related"));

  free(seqA);
  free(seqB);

  return 0;
}

