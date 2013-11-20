CC=gcc
CFLAGS=-g -O3 -ansi -Wall -Wextra -pedantic -lpthread -I.
LDLIBS=-lds -lpthread -lopt

COMPRESS_OBJS=align.o blosum62.o blosum62_matrix.o \
							coarse.o compressed.o compression.o database.o fasta.o flags.o \
							seeds.o seq.o util.o
COMPRESS_HEADERS=align.o blosum62.h coarse.h compressed.h compression.h \
								 database.h fasta.h flags.h seeds.h seq.h util.h

compress: blosum62.h cablastp-compress

cablastp-compress: cablastp-compress.o $(COMPRESS_HEADERS) $(COMPRESS_OBJS)

cablastp-compress.o: cablastp-compress.c blosum62.h database.h fasta.h

align.o: align.c align.h blosum62.h
blosum62.o: blosum62.c blosum62.h
blosum62_matrix.o: blosum62_matrix.c blosum62.h
coarse.o: coarse.c coarse.h seq.h
compressed.o: compressed.c compressed.h
compression.o: compression.c compression.h align.h coarse.h compressed.h seq.h
database.o: database.c database.h coarse.h compressed.h
fasta.o: fasta.c fasta.h util.h
flags.o: flags.c flags.h util.h
seeds.o: seeds.c seeds.h
seq.o: seq.c seq.h
uitl.o: util.c util.h

blosum62_matrix.c: ../scripts/mkBlosum
	../scripts/mkBlosum > blosum62_matrix.c

tests: test-extension test-nw test-ungapped

# test-extension: tests/test-extension.o align.o blosum62.o blosum62_matrix.o \ 
								# compression.o 
	# $(CC) $(LDFLAGS) \ 
		# tests/test-extension.o align.o blosum62.o blosum62_matrix.o compression.o \ 
		# $(LDLIBS) \ 
		# -o test-extension 

test-extension: tests/test-extension.o $(COMPRESS_OBJS)
	$(CC) $(LDFLAGS) \
		tests/test-extension.o $(COMPRESS_OBJS) \
		$(LDLIBS) \
		-o test-extension

test-nw: tests/test-nw.o align.o blosum62.o blosum62_matrix.o
	$(CC) $(LDFLAGS) \
		tests/test-nw.o align.o blosum62.o blosum62_matrix.o \
		$(LDLIBS) \
		-o test-nw

test-ungapped: tests/test-ungapped.o align.o blosum62.o blosum62_matrix.o
	$(CC) $(LDFLAGS) \
		tests/test-ungapped.o align.o blosum62.o blosum62_matrix.o \
		$(LDLIBS) \
		-o test-ungapped

clean:
	rm -f *.o tests/*.o
	rm -f cablastp-compress
	rm -f test-*

loc:
	find ./ -name '*.[ch]' -print0 | xargs -0 wc -l

