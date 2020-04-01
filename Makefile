
CFLAGS= -I. -O3
CFLAGS_DEBUG= -I. -g -O0
LIBS=-lm
CC=mpicc

evolution: evolution_mpi.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

debug: evolution_mpi.c
	$(CC) $(CFLAGS_DEBUG) -DDEBUG $< -o $@ $(LIBS)

clean:
	rm -f evolution debug

