# Makefile for vortex_id program
# contain build directives for main binaries and tests

# Compiler Choice : Gnu C Compiler ; Can change for icc
# optional CFLAGS: -Wconversion -Wall
CC = gcc
INC_DIR = -Isrc -Isrc/inih
CFLAGS = -c -Wall -O3 $(INC_DIR) 
LIBS = -lm -lgsl -lgslcblas

default: main
all: main othermain tests

othermain: mainMultiRunHistogram mainMultiRunRecursive mainMultiRunThreshold mainMultiRun2ndLamb mainMultiRunHistoShear

main: obj/main.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/mt64.o obj/vortexExtraction.o obj/ini.o obj/inputManager.o obj/essayHandler.o obj/stencilExtended.o obj/vortexExtractionExtend.o
	$(CC) -o bin/main $^ $(LIBS)

main%: obj/main%.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/mt64.o obj/vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/$@ $^ $(LIBS)

sampledVortexEssayFull: obj/sampledVortexEssayFull.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/mt64.o obj/vortexExtraction.o obj/ini.o obj/inputManager.o obj/essayHandler.o obj/stencilExtended.o obj/vortexExtractionExtend.o
	$(CC) -o bin/sampledVortexEssayFull $^ $(LIBS)

openFoamEssay: obj/openFoamEssay.o obj/lambdaInit.o obj/floodFill.o obj/vortexGen.o obj/mt64.o obj/vortexExtraction.o obj/ini.o obj/inputManager.o obj/essayHandler.o obj/stencilExtended.o obj/vortexExtractionExtend.o obj/preprocessing.o
	$(CC) -o bin/openFoamEssay $^ $(LIBS)

foamSlicer: obj/foamSlicer.o obj/ini.o obj/preprocessing.o obj/inputManager.o obj/floodFill.o obj/essayHandler.o obj/stencilExtended.o obj/lambdaInit.o obj/vortexGen.o obj/mt64.o obj/vortexExtraction.o obj/vortexExtractionExtend.o
	$(CC) -o bin/foamSlicer $^ $(LIBS)

test_lambOseenNSecGrad: obj/floodFill.o obj/lambdaInit.o obj/stencilExtended.o obj/test_lambOseenNSecGrad.o
	$(CC) -o bin/test_lambOseenNSecGrad $^ $(LIBS)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o $@

obj/test_lambOseenNSecGrad.o: src/tests/test_lambOseenNSecGrad.c
	$(CC) $(CFLAGS) src/tests/test_lambOseenNSecGrad.c -o obj/test_lambOseenNSecGrad.o

%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o obj/$@

%.o: src/tests/%.c
	$(CC) $(CFLAGS) $^ -o obj/$@

obj/mt64.o: src/mt19937-64.c src/mt64.h
	$(CC) $(CFLAGS) src/mt19937-64.c -o $@

mt64.o: src/mt19937-64.c src/mt64.h
	$(CC) $(CFLAGS) src/mt19937-64.c -o $@

obj/ini.o: src/inih/ini.c
	$(CC) $(CFLAGS) $< -o $@

ini.o: src/inih/ini.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f bin/* obj/* data/*.txt data/*.tex data/*.eps data/*.pdf data/*.aux data/*.log
