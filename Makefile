CC=gcc

all:
	$(CC) -o bin/rkseq rkseq.c
	$(CC) -fopenmp -o bin/rkpar rkpar.c

clean:
	rm -f bin/rkseq
	rm -r bin/rkpar
