P= main
OBJECTS= rbinom.o mt.o

CFLAGS= -std=gnu99 -Wall -O3

all: main

$(P): simulate.c $(OBJECTS)
	gcc $(CFLAGS) $< $(OBJECTS) -lm -o $@

%.o: %.c
	gcc $(CFLAGS) -c $< -o $@

clean:
	rm -f $(P) $(OBJECTS)
