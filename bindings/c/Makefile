# HPRMAT C Binding Makefile

CC = gcc
CFLAGS = -O3 -Wall
LDFLAGS = -L../../lib -lhprmat -lm

# Platform detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    RUN_ENV = DYLD_LIBRARY_PATH=../../lib
else
    RUN_ENV = LD_LIBRARY_PATH=../../lib
endif

.PHONY: all clean test

all: example_alpha_alpha

example_alpha_alpha: example_alpha_alpha.c hprmat.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

test: example_alpha_alpha
	$(RUN_ENV) ./example_alpha_alpha

clean:
	rm -f example_alpha_alpha
