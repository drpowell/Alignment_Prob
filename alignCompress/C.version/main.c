
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "misc.h"

#define INFINITY      1E20         // Maybe not, but close ;)

#define ALPHA_SIZE    4

char int2char(int i)
{
  static char a[ALPHA_SIZE] = "atgc";
  if (i<0 || i>ALPHA_SIZE) {
    fprintf(stderr, "Bad param int2char(%d)\n", i);
    exit(-2);
  }
  return a[i];
}

inline int char2int(char c)
{
  switch (tolower(c)) {
  case 'a':
    return 0;
  case 't':
    return 1;
  case 'g':
    return 2;
  case 'c':
    return 3;
  }
  fprintf(stderr, "Unknown char: %c\n", c);
  exit(-2);
}


// -----------------------------------------------------------------------------------
// Sequence models.
//
// res[] is really a 2d array:  double res[strlen][ALPHA_SIZE]
// that contains the prob for each character at each position along a sequence
// This must be computed for both sequences

void fit_uni(int len, char *str, double *res)
{
  int i,j;
  for (i=0; i<strlen(str); i++)
    for (j=0; j<ALPHA_SIZE; j++)
      res[i*ALPHA_SIZE + j] = log2(ALPHA_SIZE);
}

void fit_markov0(int len, char *str, double *res)
{
  int counts[ALPHA_SIZE], total;
  int i,j;
  for (i=0; i<ALPHA_SIZE; i++) counts[i]=1;

  for (i=0; i<strlen(str); i++) counts[ char2int(str[i]) ]++;
  total = ALPHA_SIZE + strlen(str);

  for (i=0; i<strlen(str); i++)
    for (j=0; j<ALPHA_SIZE; j++)
      res[i*ALPHA_SIZE + j] = -log2(1.0 * counts[j]/total);
}

void fit_markov1(int len, char *str, double *res)
{
  int counts[ALPHA_SIZE*ALPHA_SIZE], total[ALPHA_SIZE];
  int i,j;
  for (i=0; i<ALPHA_SIZE*ALPHA_SIZE; i++) counts[i]=1;
  for (i=0; i<ALPHA_SIZE; i++) total[i]=ALPHA_SIZE;

  for (i=1; i<strlen(str); i++) {
    counts[ char2int(str[i]) + ALPHA_SIZE * char2int(str[i-1]) ]++;
    total[char2int(str[i-1])]++;
  }

  for (i=0; i<strlen(str); i++)
    for (j=0; j<ALPHA_SIZE; j++) {
      double d = log2(ALPHA_SIZE);
      if (i>0)
	d = -log2(1.0 * counts[j + ALPHA_SIZE * char2int(str[i-1])]/total[char2int(str[i-1])]);
      res[i*ALPHA_SIZE + j] = d;
    }
}

// -----------------------------------------------------------------------------------


// Ok these are the FSM costs    (taken from java version)
double match_cost  = 0.5476268258420525;
double change_cost = 1.6626638258420525;

double diag_fromD  = 0.16198884271367536;
double start_fromD = 4.235036668555728;

double diag_fromI  = 0.8317598880466708;
double start_fromI = 4.904807713888723;
double cont_fromI  = 1.3048077138887235;

enum cIndex {match_i, change_i, diag_fromD_i, start_fromD_i, diag_fromI_i, start_fromI_i, cont_fromI_i, num_counts};

struct cell {
  double d,h,v;
  double counts[num_counts];
};

double cell_val(struct cell *c)
{
  return logplus(c->d, logplus(c->h, c->v));
}

void init_counts(double c[])
{
  int i;
  for (i=0; i<num_counts; i++)
    c[i] = 0.5;
}

double doAlign(char *seqA, char *seqB)
{
  int lenA = strlen(seqA);
  int lenB = strlen(seqB);

  struct cell D[2][lenB+1];
  int i,j;

  double final = INFINITY;

  double seqA_enc[lenA][ALPHA_SIZE], seqB_enc[lenB][ALPHA_SIZE];
  double seqA_cum[lenA+1], seqB_cum[lenB+1];

  //fit_uni(lenA, seqA, (double*)seqA_enc);
  //fit_uni(lenB, seqB, (double*)seqB_enc);

  //fit_markov0(lenA, seqA, (double*)seqA_enc);
  //fit_markov0(lenB, seqB, (double*)seqB_enc);

  fit_markov1(lenA, seqA, (double*)seqA_enc);
  fit_markov1(lenB, seqB, (double*)seqB_enc);
  

  // Calculate the cumulative encodings
  seqA_cum[0] = seqB_cum[0] = 0;
  for (i=0; i<lenA; i++) seqA_cum[i+1] = seqA_enc[i][char2int(seqA[i])] + seqA_cum[i];
  for (i=0; i<lenB; i++) seqB_cum[i+1] = seqB_enc[i][char2int(seqB[i])] + seqB_cum[i];


  /*  
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
  */
  



  // Clear first row
  for (j=0; j<lenB+1; j++) {
    D[0][j] = ((struct cell) {INFINITY, INFINITY, INFINITY});
    init_counts(D[0][j].counts);
  }

  for (i=0; i<lenA+1; i++) {

    // Clear next row
    for (j=0; j<lenB+1; j++) {
      D[(i+1)%2][j] = ((struct cell) {INFINITY, INFINITY, INFINITY});
      init_counts(D[(i+1)%2][j].counts);
    }
    
    for (j=0; j<lenB+1; j++) {
      
      // Add the jump _into_ D[i][j]
      D[i%2][j].d = logplus(D[i%2][j].d, seqA_cum[i]+seqB_cum[j]);

      // D[i][j] is now fully computed. Do it's jump _out_ into final cell
      final = logplus(final, 
		      cell_val(&(D[i%2][j])) + (seqA_cum[lenA]-seqA_cum[i]) + (seqB_cum[lenB]-seqB_cum[j]));

      //printf("D[%d][%d] d=%f h=%f v=%f     final=%f\n", i,j, D[i%2][j].d, D[i%2][j].h, D[i%2][j].v, final);


      // And finally outputs of D[i][j] into the 3 cells that depend on it

      // Horizontal step
      if (j<lenB) {
	double char_cost = seqB_enc[j][char2int(seqB[j])];
	double v1 = D[i%2][j].d + start_fromD + char_cost;
	double v2 = D[i%2][j].v + start_fromI + char_cost;
	double v3 = D[i%2][j].h + cont_fromI  + char_cost;
	D[i%2][j+1].h = logplus(D[i%2][j+1].h,
				logplus(v1, logplus(v2, v3)));
      }

      // Vertical step
      if (i<lenA) {
	double char_cost = seqA_enc[i][char2int(seqA[i])];
	double v1 = D[i%2][j].d + start_fromD + char_cost;
	double v2 = D[i%2][j].v + cont_fromI  + char_cost;
	double v3 = D[i%2][j].h + start_fromI + char_cost;
	D[(i+1)%2][j].v = logplus(D[(i+1)%2][j].v,
				  logplus(v1, logplus(v2, v3)));
      }

      // Diagonal step
      if (i<lenA  &&  j<lenB) {
	double diagCost, v1, v2, v3;
	if (seqA[i]==seqB[j]) {
	  // Do: P(match) * ( P(char a) + P(char b) ) / 2
	  diagCost = match_cost + logplus(seqA_enc[i][ char2int(seqA[i]) ], seqB_enc[j][ char2int(seqB[j]) ]) + 1;
	} else {
	  // Do: P(change) * P(char a) * P(char b) * 0.5 * (1/(1-P(char b)) + 1/(1-P(char a)))
	  double aN = exp2(-seqA_enc[i][ char2int(seqB[j]) ]);
	  double bN = exp2(-seqB_enc[j][ char2int(seqA[i]) ]);
	  double norm = -log2( 1/(1-aN) + 1/(1-bN) );
	  diagCost = change_cost + seqA_enc[i][ char2int(seqA[i]) ] + seqB_enc[j][ char2int(seqB[j]) ] + 1 + norm;
	}

	v1 = D[i%2][j].d + diag_fromD + diagCost;
	v2 = D[i%2][j].v + diag_fromI + diagCost;
	v3 = D[i%2][j].h + diag_fromI + diagCost;
	D[(i+1)%2][j+1].d = logplus(D[(i+1)%2][j+1].d,
				    logplus(v1, logplus(v2, v3)));
      }

    }
  }

  {
    // Add a cost for the length of the alignment.
    // Assume we know the length of the sequences. Need to encode the start and end
    // of the alignment.  Assume uniform over all positions.  Have 4 cut-points to
    // encode, but only need to encode 3 (kind of).
    // So encode 2 from the shorter sequence, and 1 from the longer.
    double l1 = (lenA < lenB ? lenA : lenB);
    double l2 = (lenA > lenB ? lenA : lenB);
    final += log2(l1) * 2 - (l1>1 ? 1 : 0);
    final += log2(l2);
  }

  printf("data len = %f\n", final);

  return seqA_cum[lenA]+seqB_cum[lenB] - final;
}




int main(int argc, char **argv)
{
  double r;
  char *seqA, *seqB;

  if (argc == 1) {
    seqA = malloc(100000);
    seqB = malloc(100000);
    fgets(seqA, 100000, stdin);
    fgets(seqB, 100000, stdin);    
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

  r = doAlign(seqA, seqB);
  printf("log odds ratio = %f\n", r);

  return 0;
}

