CC = gcc
#CFLAGS = -Wall -O3 -ffast-math
CFLAGS = -Wall -g

LIBS = -lm

all: m-align markov_pred

m-align: m-align.o markov.o common.o
	$(CC) $(CFLAGS) m-align.o markov.o common.o -o m-align $(LIBS)

markov_pred: markov_pred.o markov.o common.o
	$(CC) $(CFLAGS) markov_pred.o markov.o common.o -o markov_pred $(LIBS)

markov_pred.o: markov_pred.c markov.h common.h Makefile
m-align.o: m-align.c markov.h common.h Makefile
markov.o: markov.c markov.h common.h Makefile
common.o: common.c common.h Makefile

clean:
	rm -f m-align *.o core *.gcov *.bbg *.bb *.da gmon.out bb.out

