main: main.o lambdaInit.o floodFill.o src/lambdaInit.h src/floodFill.h
	gcc -o main obj/main.o obj/lambdaInit.o obj/floodFill.o -lm

main.o: src/main.c 
	gcc -c src/main.c -o obj/main.o

lambdaInit.o: src/lambdaInit.c src/lambdaInit.h
	gcc -c src/lambdaInit.c -o obj/lambdaInit.o

floodFill.o: src/floodFill.c src/floodFill.h
	gcc -c src/floodFill.c -o obj/floodFill.o
