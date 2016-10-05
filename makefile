# Makefile for vortex_id program
# contain build directives for main binaries and tests

# Compiler Choice : Gnu C Compiler ; Can change for icc
# optional CFLAGS: -Wconversion -Wall
CC = gcc
INC_DIR = -Isrc -Isrc/inih -Isrc/jhu -Isrc/jhtdb
CFLAGS = -c -Wall -O3 $(INC_DIR) 
LIBS = -lm -lgsl -lgslcblas

default: main
all: sampledVortexEssayFull openFoamEssay foamSlicer customVortices foamStatistics foamKinematics vortexHistogram fieldAverages fieldSpectra

othermain: mainMultiRunHistogram mainMultiRunRecursive mainMultiRunThreshold mainMultiRun2ndLamb mainMultiRunHistoShear

main: main.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

main%: main%.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o src/lambdaInit.h src/floodFill.h src/vortexGen.h src/vortexExtraction.h src/mt64.h
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

sampledVortexEssayFull: sampledVortexEssayFull.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

openFoamEssay: openFoamEssay.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

foamSlicer: foamSlicer.o ini.o preprocessing.o inputManager.o floodFill.o essayHandler.o stencilExtended.o lambdaInit.o vortexGen.o mt64.o vortexExtraction.o vortexExtractionExtend.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

customVortices: floodFill.o lambdaInit.o stencilExtended.o customVortices.o  vortexExtraction.o vortexExtractionExtend.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

foamStatistics: foamStatistics.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

foamKinematics: foamKinematics.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

foamKinematicsSmooth: foamKinematicsSmooth.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

fieldAverages: fieldAverages.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

fieldSpectra: fieldSpectra.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

vortexHistogram: vortexHistogram.o lambdaInit.o floodFill.o vortexGen.o mt64.o vortexExtraction.o ini.o inputManager.o essayHandler.o stencilExtended.o vortexExtractionExtend.o preprocessing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

jhtdbRawDownload: ini.o inputManager.o soapC.o soapClient.o stdsoap2.o turblib.o jhtdbRawDownload.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

test_XtoXbuffMirror: test_XtoXbuffMirror.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

test_uTouBuffMirror: test_uTouBuffMirror.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

test_gaussianFilterUniform: test_gaussianFilterUniform.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

test_gaussianFilterNonUniform: test_gaussianFilterNonUniform.o fieldSmoothing.o
	$(CC) -o bin/$@ $(patsubst %.o,obj/%.o,$^) $(LIBS)

%.o: src/%.c
	$(CC) $(CFLAGS) $^ -o obj/$@

%.o: src/tests/%.c
	$(CC) $(CFLAGS) $^ -o obj/$@

%.o: src/jhtdb/%.c
	$(CC) $(CFLAGS) $^ -o obj/$@

ini.o: src/inih/ini.c
	$(CC) $(CFLAGS) $< -o obj/$@

mt64.o: src/mt19937-64.c
	$(CC) $(CFLAGS) src/mt19937-64.c -o obj/$@
	
clean:
	rm -f bin/* obj/* data/*.txt data/*.tex data/*.eps data/*.pdf data/*.aux data/*.log
