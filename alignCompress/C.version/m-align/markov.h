/*  m-align
 *  Copyright (c) 2002 David Powell <david@drp.id.au>
 *
 * 
 *
    This file is part of m-align.
*/

void markov_init(int alphaSize, int order);
void markov_free();
void markov_fit(int len, unsigned char *str);
void markov_load(char *fname);
void markov_save(FILE *f);
void markov_predict(int len, unsigned char *str, DOUBLE *res);
DOUBLE markov_entropy();

