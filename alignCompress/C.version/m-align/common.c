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
#include <time.h>

#include "common.h"


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
// read_fasta - read a FASTA format file.
char *read_fasta(char *fname)
{
  char *seq;
  char buf[1024];
  int upto = 0;
  char *name = NULL;

  FILE *f = fopen(fname, "r");
  if (!f) {
    fprintf(stderr, "Unable to open sequence file '%s'.  Exiting...\n",fname);
    exit(1);
  }

  seq = malloc(MAX_SEQ_LEN);
  assert(seq && "Out of memory");

  // Find seq name first
  while (!feof(f)) {
    int i;
    if (!fgets(buf, 1024, f))
      break;
    if (buf[0] != '>') continue;
    for (i=1; buf[i] && buf[i]==' '; i++)
      ;
    if (buf[strlen(buf)-1] == '\n')
      buf[strlen(buf)-1] = 0;
    name = strdup(buf+i);
    break;
  }

  if (!name) {
    fprintf(stderr, "Unable to find fasta header line in file '%s'. Exiting...\n", fname);
    exit(1);
  }

  while (!feof(f)) {
    int len;
    int i;
    if (!fgets(buf, 1024, f))
      break;
    len = strlen(buf);
    for (i=0; i<len; i++)
      if (isalpha(buf[i])) {
        assert(upto < MAX_SEQ_LEN && "Sequence too long");
        seq[upto++] = buf[i];
      }
  }
  seq[upto] = 0;
  fprintf(stderr, "Sequence '%s': length=%d.\n", name, upto);
  //fprintf(stderr, ":%s:\n",seq);
  free(name);
  fclose(f);
  return seq;
}

