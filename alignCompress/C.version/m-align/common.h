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

#define VERSION "0.8"

#define DOUBLE double
#define DO_COUNTS


#define ALPHA_SIZE    4
const static char alphabet[ALPHA_SIZE+1] = "atgc";
//#define ALPHA_SIZE    24
//static char alphabet[ALPHA_SIZE+1] = "arndcqeghilkmfpstwyvbzx*";

#define MAX_SEQ_LEN   100000




 
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


int char2int(char c);
char int2char(int i);
void strict_DNA_seq(char *str, int len);
char *read_fasta(char *fname);



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
