all: rbinom.o

rbinom.o: rbinom.c mt.o
	gcc -std=gnu99 -Wall -O2 rbinom.c mt.o -lm -c

clean:
	rm -f rbinom.o
