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
#include <getopt.h>

#include "common.h"
#include "markov.h"

void usage(char *exeName)
{
  fprintf(stderr, "Usage: %s [-mfh] <seq>\n", exeName);
  fprintf(stderr, "  Where <seq> is a sequence file in FASTA format.\n");
  fprintf(stderr, "       -m <n>     Use a n'th order Markov Model (default 0)\n");
  fprintf(stderr, "       -f <file>  Read the Markov Model parameters from <file>\n");
  fprintf(stderr, "       -s <file>  Save the Markov Model parameters to <file>\n");
  fprintf(stderr, "       -h         This help\n");
  exit(-1);
}


int main(int argc, char **argv)
{
  char *exeName = argv[0];
  char *seqA;
  int lenA;
  int markovOrder = 0;
  char *markovFile = NULL;
  char *markovSaveFile = NULL;

  while (1) {
    int c = getopt(argc, argv, "m:f:s:h");
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
    case 'h':
    default:
      usage(exeName);
    }
  }
	
  argc -= optind-1;
  argv += optind-1;

  if (argc != 2) {
    usage(exeName);
  } else {
    seqA = read_fasta(argv[1]);
  } 

  lenA = strlen(seqA);

  printf("# Character prediction probability for FASTA file '%s'\n", argv[1]);
  printf("# Markov order = %d\n", markovOrder);
  printf("# Column order = [%s]\n", alphabet);

  {
    int i,j;
    unsigned char seqA_i[lenA];
    DOUBLE seqA_enc[lenA][ALPHA_SIZE];

    // Convert DNA sequence to only an A G C or T
    strict_DNA_seq(seqA, lenA);

    // First convert strings to numbers representing the characters
    for (i=0; i<lenA; i++) seqA_i[i] = char2int(seqA[i]);
  

    markov_init(ALPHA_SIZE, markovOrder);
    if (markovFile)
      markov_load(markovFile);
    else
      markov_fit(lenA, seqA_i);
    
    markov_predict(lenA, seqA_i, (DOUBLE*)seqA_enc);


    for (i=0; i<lenA; i++) {
      for (j=0; j<ALPHA_SIZE; j++) {
	printf("%f ", exp2(-seqA_enc[i][j]));
      }
      printf("\n");
    }

    if (markovSaveFile) {
      FILE *f = fopen(markovSaveFile, "w");
      if (!f) {
        fprintf(stderr, "Unable to open file '%s' for writing.\n", markovSaveFile);
      } else {
        fprintf(stderr, "Saving Markov Model parameters to file '%s'\n", markovSaveFile);
        markov_save(f);
      }
    }

  }

  return 0;
}
