#CC = gcc
#CFLAGS = -O3 -ffast-math -pg -fprofile-arcs -ftest-coverage -g -ax
#CFLAGS = -O3 -ffast-math

CC = gcc-3.0
#CFLAGS = -std=gnu99 -ggdb -Wall 
#CFLAGS = -std=gnu99 -O2 -ffast-math -pg -fprofile-arcs -ftest-coverage -g -ax
CFLAGS = -Wall -std=gnu99 -O3 -ffast-math
#CFLAGS = -std=gnu99 -O3
LIBS = -lm

OBJS = main.o markov.o

m-align: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o m-align $(LIBS)

main.o: main.c main.h Makefile
markov.o: markov.c main.h Makefile

clean:
	rm -f m-align *.o core *.gcov *.bbg *.bb *.da gmon.out bb.out

