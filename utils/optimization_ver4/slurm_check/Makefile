MPICC ?= mpicc
CFLAGS ?= -O2 -Wall

all: spawn_child

spawn_child: spawn_child.c
	$(MPICC) $(CFLAGS) -o $@ $<

clean:
	rm -f spawn_child

.PHONY: all clean
