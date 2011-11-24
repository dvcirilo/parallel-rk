CC=gcc

all:
	$(CC) -o bin/rkseq rkseq.c

clean:
	rm -f bin/rkseq
