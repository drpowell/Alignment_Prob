CC = gcc
CFLAGS = -ggdb -Wall
LIBS = -lm

OBJS = main.o

foo: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o foo $(LIBS)

clean:
	rm -f *.o core

