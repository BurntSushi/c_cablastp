CC=gcc
CFLAGS=-O3 -ansi -Wall -Wextra -pedantic
HEADERS=src/fasta.h
OBJS=src/main.o src/fasta.o
LDFLAGS=-I. -L. -lds

cablastp: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o cablastp

src/main.o: src/main.c
	$(CC) $(CFLAGS) -c $< -o $@

src/%.o: src/%.c src/%.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f src/*.o
	rm cablastp

