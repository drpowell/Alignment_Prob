// -----------------------------------------------------------------------------------
// Sequence models.
//
// res[] is really a 2d array:  DOUBLE res[strlen][ALPHA_SIZE]
// that contains the prob for each character at each position along a sequence
// This must be computed for both sequences

void fit_uni(int len, char *str, DOUBLE *res)
{
  int i,j;
  for (i=0; i<len; i++)
    for (j=0; j<ALPHA_SIZE; j++)
      res[i*ALPHA_SIZE + j] = log2(ALPHA_SIZE);
}

void fit_markov0(int len, char *str, DOUBLE *res)
{
  int counts[ALPHA_SIZE], total;
  int i,j;
  for (i=0; i<ALPHA_SIZE; i++) counts[i]=1;

  for (i=0; i<len; i++) counts[ (int)str[i] ]++;
  total = ALPHA_SIZE + len;

  for (i=0; i<len; i++)
    for (j=0; j<ALPHA_SIZE; j++)
      res[i*ALPHA_SIZE + j] = -log2(1.0 * counts[j]/total);
}

void fit_markov1(int len, char *str, DOUBLE *res)
{
  int counts[ALPHA_SIZE*ALPHA_SIZE], total[ALPHA_SIZE];
  int i,j;
  for (i=0; i<ALPHA_SIZE*ALPHA_SIZE; i++) counts[i]=1;
  for (i=0; i<ALPHA_SIZE; i++) total[i]=ALPHA_SIZE;

  for (i=1; i<len; i++) {
    counts[ (int)str[i] + ALPHA_SIZE * (int)str[i-1] ]++;
    total[(int)str[i-1]]++;
  }

  for (i=0; i<len; i++)
    for (j=0; j<ALPHA_SIZE; j++) {
      DOUBLE d = log2(ALPHA_SIZE);
      if (i>0)
	d = -log2(1.0 * counts[j + ALPHA_SIZE * (int)str[i-1]]/total[(int)str[i-1]]);
      res[i*ALPHA_SIZE + j] = d;
    }
}

