CC=gcc
CFLAGS=-g -O3 -ansi -Wall -Wextra -pedantic -lpthread
HEADERS=src/fasta.h src/util.h
OBJS=src/fasta.o src/main.o src/util.o
LDFLAGS=-I. -L. -lds -lpthread

cablastp: $(OBJS) src/blosum62.h
	$(CC) $(CFLAGS) $(OBJS) -o cablastp

src/main.o: src/main.c src/blosum62.h src/fasta.h
	$(CC) $(CFLAGS) -c src/main.c -o src/main.o

src/fasta.o: src/fasta.c src/fasta.h src/util.h
	$(CC) $(CFLAGS) -c src/fasta.c -o src/fasta.o

src/uitl.o: src/util.h src/util.h
	$(CC) $(CFLAGS) -c src/util.c -o src/util.o

src/blosum62.h: scripts/mkBlosum
	scripts/mkBlosum > src/blosum62.h

clean:
	rm -f src/*.o
	rm -f cablastp

