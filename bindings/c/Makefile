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

all: test_ex1
	@echo ""
	@echo "Build complete! To run:"
	@echo "  $(RUN_ENV) ./test_ex1"

test_ex1: test_ex1.c hprmat.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

test: test_ex1
	$(RUN_ENV) ./test_ex1

clean:
	rm -f test_ex1
