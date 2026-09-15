CC ?= cc
CFLAGS ?= -O2 -Wall
main: map_arch
map_arch: map_arch.v2.c
	$(CC) $(CFLAGS) -o $@ $< -lz
.PHONY: main
