CC=gcc
CFLAGS=-g -O3 -ansi -Wall -Wextra -pedantic -lpthread -I.
LDFLAGS=-L.
LDLIBS=-lds -lpthread -lopt

COMPRESS_OBJS=cablastp-compress.o fasta.o util.o
COMPRESS_HEADERS=blosum62.h fasta.h util.h

compress: blosum62.h cablastp-compress

cablastp-compress: $(COMPRESS_HEADERS) $(COMPRESS_OBJS)

main.o: main.c blosum62.h fasta.h

fasta.o: fasta.c fasta.h util.h

uitl.o: util.c util.h

blosum62.h: ../scripts/mkBlosum
	../scripts/mkBlosum > blosum62.h

clean:
	rm -f *.o
	rm -f cablastp-compress

