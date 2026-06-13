MPICC ?= mpicc
CFLAGS ?= -O2 -Wall

all: worker

worker: worker.c
	$(MPICC) $(CFLAGS) -o $@ $<

clean:
	rm -f worker
	rm -rf out

.PHONY: all clean
