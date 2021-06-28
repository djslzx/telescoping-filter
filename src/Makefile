#compiler
CC=gcc	
#CC = clang	

#put any desired compilation flags here.  Feel free to remove O2 or change to O3 or Ofast	
#make sure to run "make clean" if this is changed	
CFLAGS=-O3 -lm	

#flags to use when running make debug	
#replaces CFLAGS	
DEBUGFLAGS=-g -O0 -lm	

#flags to use when running gprof
PROFFLAGS=-pg -O0 -lm

#operating system (for Max)
OS := $(shell uname)
ifeq ($(OS), Darwin)
CFLAGS += -Wno-unused-command-line-argument -Wformat=0
PROFFLAGS += -Wno-unused-command-line-argument -Wformat=0
else
endif

#files to be kept up-to-date	
OBJ = test.o filter.o murmur3.o set.o remote.o arithmetic_coding.o dict.o

#include any files beyond .c files that depend on .o object files here	
#(this likely only includes headers)	
DEPS = filter.h murmur3.h set.h remote.h arithmetic_coding.h dict.h remainder.h constants.h

#only need test.out to build 'all' of project
all: test.out	

unit_test: test.out
	./test.out 1

integration_test: test.out
	./test.out 2

test_all: test.out
	./test.out 0 && ./test.out "../data/chicagoA64.txt" 605006 100 0.95

#runs tests (Chicago A)
#exe file num-unique-lines A/S load
test: test.out
	./test.out "../data/chicagoA64.txt" 605006 100 0.95

#runs tests with memory leak check
vgtest: test.out
	valgrind --leak-check=full ./test.out 2

#adds flags (set above) for make debug	
#make sure to run "make clean" if no changes to source files	
debug: CFLAGS=$(DEBUGFLAGS)
debug: test.out

#profile with gprof
prof: CFLAGS=$(PROFFLAGS)
prof: test.out

#this rule says that every .o file needs to be compiled using the corresponding .c file	
%.o: %.c $(DEPS)	
	$(CC) -c -o $@ $< $(CFLAGS)	

#this rule links the object files together and stores the output in test.out	
test.out: $(OBJ)	
	$(CC) -o $@ $^ $(CFLAGS)	

#a possibly-sloppy way to undo making: remove all object files	
clean: 	
	rm $(OBJ)
