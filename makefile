# Makefile for vortex_id program
# contain build directives for main binaries and tests

# Compiler Choice : Gnu C Compiler ; Can change for icc
CC = gcc
INC_DIR = src
CFLAGS = -c -O3 -I$(INC_DIR)
LIBS = -lm -lgsl -lgslcblas

default: main

all: main test_floodFill test_lambOseenShear0 test_lambOseenShear1 test_lambOseenShear2 test_lambOseenShear3 test_addSingleOseen test_genLOseenUniformList test_genLOseenBinaryList test_vortexQuickSort test_vortexExtraction0 test_vortexExtraction1 test_vortexExtraction2 test_vortexExtraction3 test_vortexSingleRun test_vortexSingleRunTime test_vortexShearSingleRunTime test_vortexMultiRun test_vortexMultiRunHistogram test_vortexExtSimple test_vortexExtRecursive

main: main.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/main obj/main.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS)

test_floodFill: floodFill.o test_floodFill.o
	$(CC) -o bin/test_floodFill obj/test_floodFill.o obj/floodFill.o $(LIBS)

test_lambOseenShear0: floodFill.o test_lambOseenShear0.o lambdaInit.o
	$(CC) -o bin/test_lambOseenShear0 obj/test_lambOseenShear0.o obj/floodFill.o obj/lambdaInit.o $(LIBS)

test_lambOseenShear1: floodFill.o test_lambOseenShear1.o lambdaInit.o
	$(CC) -o bin/test_lambOseenShear1 obj/test_lambOseenShear1.o obj/floodFill.o obj/lambdaInit.o $(LIBS)

test_lambOseenShear2: floodFill.o test_lambOseenShear2.o lambdaInit.o
	$(CC) -o bin/test_lambOseenShear2 obj/test_lambOseenShear2.o obj/floodFill.o obj/lambdaInit.o $(LIBS)

test_lambOseenShear3: floodFill.o test_lambOseenShear3.o lambdaInit.o
	$(CC) -o bin/test_lambOseenShear3 obj/test_lambOseenShear3.o obj/floodFill.o obj/lambdaInit.o $(LIBS)

test_addSingleOseen: test_addSingleOseen.o lambdaInit.o floodFill.o
	$(CC) -o bin/test_addSingleOseen obj/test_addSingleOseen.o obj/floodFill.o obj/lambdaInit.o $(LIBS)

test_genLOseenUniformList: test_genLOseenUniformList.o mt64.o vortexGen.o
	$(CC) -o bin/test_genLOseenUniformList obj/test_genLOseenUniformList.o obj/mt64.o obj/vortexGen.o $(LIBS)

test_genLOseenBinaryList: test_genLOseenBinaryList.o mt64.o vortexGen.o
	$(CC) -o bin/test_genLOseenBinaryList obj/test_genLOseenBinaryList.o obj/mt64.o obj/vortexGen.o $(LIBS)

test_vortexQuickSort: test_vortexQuickSort.o mt64.o vortexGen.o vortexExtraction.o lambdaInit.o floodFill.o
	$(CC) -o bin/test_vortexQuickSort obj/test_vortexQuickSort.o obj/mt64.o obj/vortexGen.o obj/vortexExtraction.o obj/lambdaInit.o obj/floodFill.o $(LIBS)

test_vortexExtraction0: test_vortexExtraction0.o floodFill.o lambdaInit.o vortexExtraction.o
	$(CC) -o bin/test_vortexExtraction0 obj/test_vortexExtraction0.o obj/floodFill.o obj/lambdaInit.o obj/vortexExtraction.o $(LIBS)

test_vortexExtraction1: test_vortexExtraction1.o floodFill.o lambdaInit.o vortexExtraction.o
	$(CC) -o bin/test_vortexExtraction1 obj/test_vortexExtraction1.o obj/floodFill.o obj/lambdaInit.o obj/vortexExtraction.o $(LIBS)

test_vortexExtraction2: test_vortexExtraction2.o floodFill.o lambdaInit.o vortexExtraction.o
	$(CC) -o bin/test_vortexExtraction2 obj/test_vortexExtraction2.o obj/floodFill.o obj/lambdaInit.o obj/vortexExtraction.o $(LIBS)

test_vortexExtraction3: test_vortexExtraction3.o floodFill.o lambdaInit.o vortexExtraction.o
	$(CC) -o bin/test_vortexExtraction3 obj/test_vortexExtraction3.o obj/floodFill.o obj/lambdaInit.o obj/vortexExtraction.o $(LIBS)

test_vortexExtraction4: test_vortexExtraction4.o floodFill.o lambdaInit.o vortexExtraction.o
	$(CC) -o bin/test_vortexExtraction4 obj/test_vortexExtraction4.o obj/floodFill.o obj/lambdaInit.o obj/vortexExtraction.o $(LIBS)

test_vortexSingleRun: test_vortexSingleRun.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexSingleRun obj/test_vortexSingleRun.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS)

test_vortexSingleRunTime: test_vortexSingleRunTime.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexSingleRunTime obj/test_vortexSingleRunTime.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

test_vortexShearSingleRunTime: test_vortexShearSingleRunTime.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexShearSingleRunTime obj/test_vortexShearSingleRunTime.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

test_vortexMultiRun: test_vortexMultiRun.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexMultiRun obj/test_vortexMultiRun.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

test_vortexMultiRunHistogram: test_vortexMultiRunHistogram.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexMultiRunHistogram obj/test_vortexMultiRunHistogram.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

test_vortexExtSimple: test_vortexExtSimple.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexExtSimple obj/test_vortexExtSimple.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

test_vortexExtRecursive: test_vortexExtRecursive.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/test_vortexExtRecursive obj/test_vortexExtRecursive.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/vortexExtraction.o obj/mt64.o $(LIBS) 

main.o: src/main.c 
	$(CC) $(CFLAGS) src/main.c -o obj/main.o

lambdaInit.o: src/lambdaInit.c src/lambdaInit.h
	$(CC) $(CFLAGS) src/lambdaInit.c -o obj/lambdaInit.o

floodFill.o: src/floodFill.c src/floodFill.h
	$(CC) $(CFLAGS) src/floodFill.c -o obj/floodFill.o

vortexGen.o: src/vortexGen.c
	$(CC) $(CFLAGS) src/vortexGen.c -o obj/vortexGen.o

vortexExtraction.o: src/vortexExtraction.c
	$(CC) $(CFLAGS) src/vortexExtraction.c -o obj/vortexExtraction.o

mt64.o: src/mt19937-64.c src/mt64.h
	$(CC) $(CFLAGS) src/mt19937-64.c -o obj/mt64.o

test_floodFill.o: src/tests/test_floodFill.c src/floodFill.h
	$(CC) $(CFLAGS) src/tests/test_floodFill.c -o obj/test_floodFill.o

test_lambOseenShear0.o: src/tests/test_lambOseenShear0.c src/floodFill.h src/lambdaInit.h
	$(CC) $(CFLAGS) src/tests/test_lambOseenShear0.c -o obj/test_lambOseenShear0.o

test_lambOseenShear1.o: src/tests/test_lambOseenShear1.c src/floodFill.h src/lambdaInit.h
	$(CC) $(CFLAGS) src/tests/test_lambOseenShear1.c -o obj/test_lambOseenShear1.o

test_lambOseenShear2.o: src/tests/test_lambOseenShear2.c src/floodFill.h src/lambdaInit.h
	$(CC) $(CFLAGS) src/tests/test_lambOseenShear2.c -o obj/test_lambOseenShear2.o

test_lambOseenShear3.o: src/tests/test_lambOseenShear3.c src/floodFill.h src/lambdaInit.h
	$(CC) $(CFLAGS) src/tests/test_lambOseenShear3.c -o obj/test_lambOseenShear3.o

test_addSingleOseen.o: src/tests/test_addSingleOseen.c src/floodFill.h src/lambdaInit.h
	$(CC) $(CFLAGS) src/tests/test_addSingleOseen.c -o obj/test_addSingleOseen.o

test_genLOseenUniformList.o: src/tests/test_genLOseenUniformList.c 
	$(CC) $(CFLAGS) src/tests/test_genLOseenUniformList.c -o obj/test_genLOseenUniformList.o

test_genLOseenBinaryList.o: src/tests/test_genLOseenBinaryList.c 
	$(CC) $(CFLAGS) src/tests/test_genLOseenBinaryList.c -o obj/test_genLOseenBinaryList.o

test_vortexQuickSort.o: src/tests/test_vortexQuickSort.c 
	$(CC) $(CFLAGS) src/tests/test_vortexQuickSort.c -o obj/test_vortexQuickSort.o

test_vortexExtraction0.o: src/tests/test_vortexExtraction0.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtraction0.c -o obj/test_vortexExtraction0.o

test_vortexExtraction1.o: src/tests/test_vortexExtraction1.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtraction1.c -o obj/test_vortexExtraction1.o

test_vortexExtraction2.o: src/tests/test_vortexExtraction2.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtraction2.c -o obj/test_vortexExtraction2.o

test_vortexExtraction3.o: src/tests/test_vortexExtraction3.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtraction3.c -o obj/test_vortexExtraction3.o

test_vortexExtraction4.o: src/tests/test_vortexExtraction4.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtraction4.c -o obj/test_vortexExtraction4.o

test_vortexSingleRun.o: src/tests/test_vortexSingleRun.c
	$(CC) $(CFLAGS) src/tests/test_vortexSingleRun.c -o obj/test_vortexSingleRun.o

test_vortexSingleRunTime.o: src/tests/test_vortexSingleRunTime.c
	$(CC) $(CFLAGS) src/tests/test_vortexSingleRunTime.c -o obj/test_vortexSingleRunTime.o

test_vortexShearSingleRunTime.o: src/tests/test_vortexShearSingleRunTime.c
	$(CC) $(CFLAGS) src/tests/test_vortexShearSingleRunTime.c -o obj/test_vortexShearSingleRunTime.o

test_vortexMultiRun.o: src/tests/test_vortexMultiRun.c
	$(CC) $(CFLAGS) src/tests/test_vortexMultiRun.c -o obj/test_vortexMultiRun.o

test_vortexMultiRunHistogram.o: src/tests/test_vortexMultiRunHistogram.c
	$(CC) $(CFLAGS) src/tests/test_vortexMultiRunHistogram.c -o obj/test_vortexMultiRunHistogram.o

test_vortexExtSimple.o: src/tests/test_vortexExtSimple.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtSimple.c -o obj/test_vortexExtSimple.o

test_vortexExtRecursive.o: src/tests/test_vortexExtRecursive.c
	$(CC) $(CFLAGS) src/tests/test_vortexExtRecursive.c -o obj/test_vortexExtRecursive.o

clean:
	rm bin/* obj/* data/*.txt data/*.tex data/*.eps
