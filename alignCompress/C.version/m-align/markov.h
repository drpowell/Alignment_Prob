/*  m-align
 *  Copyright (c) 2002 David Powell
 *
 * 
 *
    This file is part of m-align.

    m-align is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    m-align is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

*/

void markov_init(int alphaSize, int order);
void markov_free();
void markov_fit(int len, unsigned char *str);
void markov_load(char *fname);
void markov_save(FILE *f);
void markov_predict(int len, unsigned char *str, DOUBLE *res);
DOUBLE markov_entropy();

