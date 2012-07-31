CC=gcc
CFLAGS=-g -O3 -ansi -Wall -Wextra -pedantic -lpthread -I.
LDFLAGS=-L.
LDLIBS=-lds -lpthread -lopt

COMPRESS_OBJS=cablastp-compress.o coarse.o compressed.o database.o fasta.o flags.o seq.o util.o
COMPRESS_HEADERS=blosum62.h coarse.h compressed.h database.h fasta.h flags.h seq.h util.h

compress: blosum62.h cablastp-compress

cablastp-compress: $(COMPRESS_HEADERS) $(COMPRESS_OBJS)

cablastp-compress.o: cablastp-compress.c blosum62.h database.h fasta.h

coarse.o: coarse.c coarse.h seq.h
compressed.o: compressed.c compressed.h
database.o: database.c database.h coarse.h compressed.h
fasta.o: fasta.c fasta.h util.h
flags.o: flags.c flags.h util.h
seq.o: seq.c seq.h
uitl.o: util.c util.h

blosum62.h: ../scripts/mkBlosum
	../scripts/mkBlosum > blosum62.h

clean:
	rm -f *.o
	rm -f cablastp-compress

