#compiler
CC=gcc
LD=ld
#CC = clang	

CFLAGS += -lm
CFLAGS += -fsanitize=undefined
#CFLAGS += -fsanitize=undefined-abort

#put any desired compilation flags here.  Feel free to remove O2 or change to O3 or Ofast	
#make sure to run "make clean" if this is changed	
RELFLAGS=$(CFLAGS) -O3

#flags to use when running make debug	
#replaces RELFLAGS	
DEBUGFLAGS=$(CFLAGS) -g -O0 -rdynamic -Wall

#flags to use when running gprof
PROFFLAGS=$(CFLAGS) -pg -O0

#operating system (for Max)
OS := $(shell uname)
ifeq ($(OS), Darwin)
RELFLAGS += -Wno-unused-command-line-argument -Wformat=0
PROFFLAGS += -Wno-unused-command-line-argument -Wformat=0
else
endif

DEPS = arcd.h constants.h macros.h murmur3.h bit_util.h remainder.h rsqf.h set.h
OBJ = arcd.o exaf.o murmur3.o bit_util.o rsqf.o set.o
ALGO = rsqf exaf utaf taf arcd

#only need test.out to build 'all' of project
all: test.out $(ALGO)

#adds flags (set above) for make debug
#make sure to run "make clean" if no changes to source files	
debug: RELFLAGS=$(DEBUGFLAGS)
debug: test.out

#profile with gprof
prof: RELFLAGS=$(PROFFLAGS)
prof: test.out

.PHONY: all clean

rsqf: rsqf.c
	$(CC) -D TEST_RSQF=1 -o rsqf rsqf.c murmur3.c bit_util.c set.c $(DEBUGFLAGS)

exaf: exaf.c
	$(CC) -D TEST_EXAF=1 -o exaf exaf.c arcd.c murmur3.c bit_util.c set.c $(DEBUGFLAGS)

utaf: utaf.c
	$(CC) -D TEST_UTAF=1 -o utaf utaf.c arcd.c murmur3.c bit_util.c set.c $(DEBUGFLAGS)

taf: taf.c
	$(CC) -D TEST_TAF=1 -o taf taf.c arcd.c murmur3.c bit_util.c set.c $(DEBUGFLAGS)

arcd: arcd.c
	$(CC) -D TEST_ARCD=1 -o arcd arcd.c $(DEBUGFLAGS)

# $@ = target name
# $^ = all prereqs

# this rule says that every .o file needs to be compiled using the corresponding .c file
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(RELFLAGS)

#this rule links the object files together and stores the output in test.out	
test.out: $(OBJ)	
	$(LD) -relocatable -o $@ $^

clean: 	
	rm test.out $(OBJ) $(ALGO)
